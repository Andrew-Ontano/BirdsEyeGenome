#!/usr/bin/env python3
"""Birds Eye Genome (BEG) core function caller and CLI entry point."""

import argparse
import sys
from collections.abc import Sequence

from src.commands.count import run_count
from src.commands.plot import run_plot


def build_parser() -> argparse.ArgumentParser:
    """Build and return the main ArgumentParser with subcommands.

    Returns:
        argparse.ArgumentParser: Configured argument parser.
    """
    parser = argparse.ArgumentParser(
        prog="beg",
        description="Birds Eye Genome: Genomic sequence signature at a distance.",
    )
    subparsers = parser.add_subparsers(
        dest="command",
        title="subcommands",
        description="Valid subcommands for Birds Eye Genome workflow",
        help="Sub-command help",
    )

    # Subparser for count
    parser_count = subparsers.add_parser(
        "count",
        help="Count k-mers across genomic sequences via sliding windows.",
    )
    parser_count.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Path to input FASTA file.",
    )
    parser_count.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        default=None,
        help="Path to output BEG file.",
    )
    parser_count.set_defaults(func=run_count)

    # Subparser for plot
    parser_plot = subparsers.add_parser(
        "plot",
        help="Plot k-mer distributions and sequence signatures.",
    )
    parser_plot.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Path to input BEG file or sequence data.",
    )
    parser_plot.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        default=None,
        help="Path to output image/figure.",
    )
    parser_plot.set_defaults(func=run_plot)

    return parser


def main(argv: Sequence[str] | None = None) -> None:
    """Main CLI entry point.

    Args:
        argv: Optional list of command line arguments. If None, sys.argv[1:] is used.
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        sys.exit(0)

    args.func(args)


if __name__ == "__main__":
    main()
