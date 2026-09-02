import argparse
from pathlib import Path
import subprocess

from diffloop.counts import run_counting


def run_deseq2(
    rawcounts,
    normalization,
    treatment,
    control,
    treatnum,
    ctrlnum,
    outdir,
    prefix,
    log2fc,
    padj
):

    rscript = (
        Path(__file__).parent /
        "r" /
        "deseq2.R"
    )

    cmd = [
        "Rscript",
        str(rscript),
        "--Pathid", str(outdir),
        "--Rawcounts", str(rawcounts),
        "--Normalization", str(normalization),
        "--Treatid", treatment,
        "--Ctrlid", control,
        "--Treatnum", str(treatnum),
        "--Ctrlnum", str(ctrlnum),
        "--outfile", prefix,
        "--log2FC", str(log2fc),
        "--padj", str(padj)
    ]

    subprocess.run(cmd, check=True)


def run_pipeline(args):

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("================================")
    print(" DiffLoop")
    print(" Differential chromatin loops")
    print("================================")

    print("\n[1/2] Extracting loop counts")

    run_counting(
        samples=args.samples,
        loopfile=args.loops,
        outdir=outdir,
        prefix=args.prefix
    )

    rawcounts = outdir / f"{args.prefix}_rawcounts.txt"

    normalization = (
        outdir /
        f"{args.prefix}_normalizationFactor.txt"
    )

    print("\n[2/2] Running DESeq2")

    run_deseq2(
        rawcounts=rawcounts,
        normalization=normalization,
        treatment=args.treatment,
        control=args.control,
        treatnum=args.treatnum,
        ctrlnum=args.ctrlnum,
        outdir=outdir,
        prefix=args.prefix,
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

    subparsers = parser.add_subparsers(
        dest="command",
        required=True
    )

    run = subparsers.add_parser(
        "run",
        help="Run complete differential loop analysis"
    )

    run.add_argument(
        "--samples",
        required=True
    )

    run.add_argument(
        "--loops",
        required=True
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
        "--outdir",
        required=True
    )

    run.add_argument(
        "--prefix",
        required=True
    )

    run.add_argument(
        "--log2fc",
        type=float,
        default=1
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
