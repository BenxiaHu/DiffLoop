import argparse
from pathlib import Path
import subprocess

from diffloop.count import get_count
from diffloop._version import __version__


def run_deseq2(
    rawcounts,
    normalization,
    treatment,
    control,
    treatnum,
    ctrlnum,
    outpath,
    outfile,
    log2fc,
    padj
):
    rscript = Path(__file__).parent / "r" / "deseq2.R"

    cmd = [
        "Rscript",
        str(rscript),
        "--Pathid", str(outpath),
        "--Rawcounts", str(rawcounts),
        "--Normalization", str(normalization),
        "--Treatid", treatment,
        "--Ctrlid", control,
        "--Treatnum", str(treatnum),
        "--Ctrlnum", str(ctrlnum),
        "--outfile", outfile,
        "--log2FC", str(log2fc),
        "--padj", str(padj)
    ]

    subprocess.run(cmd, check=True)


def run_pipeline(args):

    outpath = Path(args.outpath)
    outpath.mkdir(parents=True, exist_ok=True)

    print("================================")
    print(" DiffLoop")
    print(" Differential chromatin loops")
    print("================================")

    print("\n[1/2] Extracting loop counts")

    rawcounts, normalization = get_count(
        inputfile=args.inputfile,
        loopfile=args.loopfile,
        outpath=outpath,
        outfile=args.outfile
    )

    print(f"Raw counts: {rawcounts}")
    print(f"Normalization factors: {normalization}")

    print("\n[2/2] Running DESeq2")

    run_deseq2(
        rawcounts=rawcounts,
        normalization=normalization,
        treatment=args.treatment,
        control=args.control,
        treatnum=args.treatnum,
        ctrlnum=args.ctrlnum,
        outpath=outpath,
        outfile=args.outfile,
        log2fc=args.log2fc,
        padj=args.padj
    )

    print("\nDiffLoop completed successfully.")


def main():

    parser = argparse.ArgumentParser(
        prog="DiffLoop",
        description=(
            "Identify differential chromatin loops "
            "between two conditions."
        )
    )

    parser.add_argument(
        "-V",
        "--version",
        action="version",
        version=f"DiffLoop {__version__}",
        help="Print version and exit"
    )

    subparsers = parser.add_subparsers(
        dest="command",
        required=True
    )

    run = subparsers.add_parser(
        "run",
        help="Run complete differential loop analysis"
    )

    run.add_argument(
        "-i",
        "--inputfile",
        dest="inputfile",
        required=True,
        help=(
            "Sample metadata TSV with columns: "
            "sample, mcool_path, resolution, expected_path"
        )
    )

    run.add_argument(
        "-l",
        "--loopfile",
        dest="loopfile",
        required=True,
        help="Loop file"
    )

    run.add_argument(
        "--treatment",
        required=True
    )

    run.add_argument(
        "--control",
        required=True
    )

    run.add_argument(
        "--treatnum",
        type=int,
        required=True
    )

    run.add_argument(
        "--ctrlnum",
        type=int,
        required=True
    )

    run.add_argument(
        "-O",
        "--outpath",
        dest="outpath",
        required=True,
        help="Output directory"
    )

    run.add_argument(
        "-o",
        "--outfile",
        dest="outfile",
        required=True,
        help="Output prefix"
    )

    run.add_argument(
        "--log2fc",
        type=float,
        default=1.0
    )

    run.add_argument(
        "--padj",
        type=float,
        default=0.05
    )

    run.set_defaults(func=run_pipeline)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
