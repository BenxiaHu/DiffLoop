#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import cooler
import os
import sys, getopt
import os.path

dir = os.path.dirname(__file__)
version_py = os.path.join(dir, "_version.py")
exec(open(version_py).read())

# ---------------------------- helpers ----------------------------

def canon_pair(row):
    """Canonicalize a loop by sorting anchors; returns (c1,s1,e1,c2,s2,e2)."""
    a = (row["chrom1"], int(row["start1"]), int(row["end1"]))
    b = (row["chrom2"], int(row["start2"]), int(row["end2"]))
    A, B = sorted([a, b], key=lambda x: (x[0], x[1], x[2]))
    return A + B

def _make_loop_id(df):
    return (df["chrom1"] + ":" + df["start1"].astype(str) + "-" + df["end1"].astype(str) +
            "|" + df["chrom2"] + ":" + df["start2"].astype(str) + "-" + df["end2"].astype(str))

def load_loops(loopfile):
    """
    Load loops, add swapped orientation, and derive a canonical intra-chrom loop index (meta).
    Returns (loops_all, meta) where:
      loops_all: all loops with swapped orientation
      meta: canonical intra-chrom loops with loop_id and distance
    """
    loops = pd.read_csv(loopfile, sep="\t", header=None)
    loops = loops.iloc[:, :6].copy()
    loops.columns = ['chrom1', 'start1', 'end1','chrom2', 'start2', 'end2']

    # add swapped orientation so both (A,B) and (B,A) match pixels
    sw = loops.rename(columns={
        'chrom1': 'chrom2', 'start1': 'start2', 'end1': 'end2',
        'chrom2': 'chrom1', 'start2': 'start1', 'end2': 'end1'
    })
    loops_all = (pd.concat([loops, sw], ignore_index=True).drop_duplicates().reset_index(drop=True))

    can = loops_all.apply(canon_pair, axis=1, result_type='expand')
    can.columns = ["c1", "s1", "e1", "c2", "s2", "e2"]
    meta = can.drop_duplicates().rename(columns={
        'c1': 'chrom1', 's1': 'start1', 'e1': 'end1',
        'c2': 'chrom2', 's2': 'start2', 'e2': 'end2'
    })

    # intra-chrom only; compute distance and loop_id
    meta = meta[meta['chrom1'] == meta['chrom2']].copy()
    meta['dist_bp'] = (meta['start2'] - meta['start1']).abs().astype(int)
    meta['loop_id'] = _make_loop_id(meta)
    meta = meta.set_index('loop_id')
    return loops_all, meta

def get_counts_for_sample(mcool_path, resolution, loops_all, meta):
    """Exact pixel counts (sum over orientations) collapsed to canonical loop_id."""
    c = cooler.Cooler(f"{mcool_path}::/resolutions/{int(resolution)}")
    # NOTE: This loads all pixels at the given resolution. Consider per-chr fetching if memory is tight.
    pix = c.matrix(balance=False, as_pixels=True, join=True)[:]
    pix = pix[['chrom1','start1','end1','chrom2','start2','end2','count']].copy()
    pix['count'] = pix['count'].fillna(0)

    m = pix.merge(loops_all, on=['chrom1','start1','end1','chrom2','start2','end2'], how='inner')

    can = m.apply(canon_pair, axis=1, result_type='expand')
    can.columns = ["c1","s1","e1","c2","s2","e2"]
    m_can = pd.concat([m, can], axis=1)
    m_can['loop_id'] = (
        m_can['c1'] + ":" + m_can['s1'].astype(str) + "-" + m_can['e1'].astype(str) +
        "|" + m_can['c2'] + ":" + m_can['s2'].astype(str) + "-" + m_can['e2'].astype(str)
    )

    counts = m_can.groupby('loop_id', as_index=True)['count'].sum()
    counts = counts.reindex(meta.index).fillna(0).astype(int)  # align to canonical set; keep zeroes
    return counts

def get_expected_for_sample_exact(expected_path, meta):
    """
    Accepts expected as TSV with either:
      (region1, region2, dist_bp, count.avg) or (chrom1, chrom2, dist_bp, count.avg)
    Filters to intra-chrom rows and returns a vector aligned to meta.index (loop_id).
    """
    e = pd.read_csv(expected_path, sep="\t", header=0)
    if {'region1','region2','dist_bp','count.avg'}.issubset(e.columns):
        e = e[e['region1'] == e['region2']].copy()
        e = e.rename(columns={'region1':'chrom1','region2':'chrom2','count.avg':'expected'})
    elif {'chrom1','chrom2','dist_bp','count.avg'}.issubset(e.columns):
        e = e[e['chrom1'] == e['chrom2']].copy()
        e = e.rename(columns={'count.avg':'expected'})
    else:
        raise ValueError(
            "expected file must contain either (region1, region2, dist_bp, count.avg) "
            "or (chrom1, chrom2, dist_bp, count.avg)"
        )

    e = e[['chrom1','chrom2','dist_bp','expected']].copy()
    e = e[e['dist_bp'] > 0]
    e['dist_bp'] = e['dist_bp'].astype(int)
    e['expected'] = e['expected'].fillna(0)

    meta_reset = meta.reset_index()[['loop_id','chrom1','start1','end1','chrom2','start2','end2','dist_bp']]
    # exact INNER merge (your rule)
    joined = meta_reset.merge(e, on=['chrom1','chrom2','dist_bp'], how='inner').set_index('loop_id')
    # re-align to full meta to guarantee matching matrices; missing expected stays NaN (will yield nf=0.0)
    expected_vec = joined['expected'].reindex(meta.index)
    return expected_vec

def read_manifest(samples_path):
    """
    Reads a TSV with columns: sample, mcool_path, resolution, expected_path
    """
    man = pd.read_csv(samples_path, sep="\t", header=0)
    need = {"sample","mcool_path","resolution","expected_path"}
    missing = need - set(man.columns)
    if missing:
        raise ValueError(f"Missing column(s) in samples TSV: {', '.join(sorted(missing))}")
    man['sample'] = man['sample'].astype(str)
    man['resolution'] = man['resolution'].astype(int)
    return man

def build_matrices(man, loops_all, meta):
    """Return (counts_df, nf_df) with identical index (loop_id) and columns (samples)."""
    counts_cols, nf_cols = {}, {}
    for _, r in man.iterrows():
        s, mc, res, ex = r['sample'], r['mcool_path'], r['resolution'], r['expected_path']
        print(f"[{s}] counting exact loop pixels…")
        counts_cols[s] = get_counts_for_sample(mc, res, loops_all, meta)

        print(f"[{s}] merging expected (exact chrom1,chrom2,dist_bp)…")
        exp_vec = get_expected_for_sample_exact(ex, meta)

        # nf = 1 / expected (your pattern). expected==0 or NaN -> nf = 0.0
        nf = 1.0 / exp_vec.replace(0, np.nan)
        nf = nf.fillna(0.0)
        nf_cols[s] = nf

    samples = man['sample'].tolist()
    counts_df = pd.DataFrame(counts_cols, index=meta.index)[samples]
    nf_df     = pd.DataFrame(nf_cols,     index=meta.index)[samples]
    return counts_df, nf_df

def write_matrices(counts_df, nf_df, out_dir, out_prefix):
    """
    Write two matched files into out_dir with prefix out_prefix:
      <out_prefix>_rawcounts.txt
      <out_prefix>_normalizationFactor.txt

    Returns
    -------
    tuple[Path, Path]
        Paths to the raw-count and normalization-factor files.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    counts_path = out_dir / f"{out_prefix}_rawcounts.txt"
    nf_path = out_dir / f"{out_prefix}_normalizationFactor.txt"

    counts_out = counts_df.copy()
    counts_out.insert(0, "loop_id", counts_out.index)
    counts_path.write_text(
        counts_out.to_csv(sep="\t", index=False, header=True)
    )

    nf_out = nf_df.copy()
    nf_out.insert(0, "loop_id", nf_out.index)
    nf_path.write_text(
        nf_out.to_csv(sep="\t", index=False, header=True)
    )

    return counts_path, nf_path


def get_count(inputfile,loopfile, outpath, outfile):
    loops_all, meta = load_loops(loopfile)
    man = read_manifest(inputfile)
    counts_df, nf_df = build_matrices(man, loops_all, meta)
    return write_matrices(counts_df, nf_df, outpath, outfile)

# ---------------------------- CLI ----------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Build loop raw count and normalization factor matrices from mcool and expected files."
    )
    parser.add_argument('-i', '--inputfile', dest='inputfile', required=True,
                        help='Sample filename (tab-separated with columns: sample, mcool_path, resolution, expected_path)')
    parser.add_argument('-l', '--loopfile', dest='loopfile', required=True,
                        help='Loop file with BIN1_*/BIN2_* columns')
    parser.add_argument('-O', '--outpath', dest='outpath', required=True,
                        help='Output directory')
    parser.add_argument('-o', '--outfile', dest='outfile', required=True,
                        help='Output prefix (no extension) for the two result files')
    parser.add_argument("-V", "--version", action="version",version="DiffLoop {}".format(__version__)\
                      ,help="Print version and exit")
    args = parser.parse_args()
    print('###Parameters:')
    print(args)
    print('###Parameters')

    get_count(args.inputfile, args.loopfile, args.outpath, args.outfile)

if __name__ == "__main__":
    main()
