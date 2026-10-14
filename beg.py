#!/usr/bin/env python3
"""Birds Eye Genome (BEG) core function caller and CLI entry point."""

import argparse
import logging
import sys

from src.commands.count import run_count
from src.commands.plot import run_plot


def configure_logging(verbose: bool) -> None:
    """Configure standard logging level and format based on verbosity flag.

    Args:
        verbose (bool): If True, log level is set to INFO; otherwise ERROR.
    """
    log_level = logging.INFO if verbose else logging.ERROR
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        force=True,
    )


def build_parser() -> argparse.ArgumentParser:
    """Build and return the main ArgumentParser with subcommands.

    Returns:
        argparse.ArgumentParser: Configured argument parser.
    """
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="Enable verbose output (sets logging level to INFO).",
    )

    parser = argparse.ArgumentParser(
        prog="beg",
        description="Birds Eye Genome: Genomic sequence signature at a distance.",
        parents=[parent_parser],
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
        parents=[parent_parser],
    )
    parser_count.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Path to input FASTA file.",
    )
    parser_count.add_argument(
        "-f",
        "--input-format",
        type=str,
        required=False,
        default="fasta",
        help="Input sequence file format (supported by BioPython).",
    )
    parser_count.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        default=None,
        help="Path to output BEG file.",
    )
    parser_count.add_argument(
        "-p",
        "--output-format",
        type=str,
        required=False,
        default="beg",
        choices=["beg", "db"],
        help="Output format. Options: tab-separated values (beg); sqlite database (db).",
    )
    parser_count.add_argument(
        "-w",
        "--window-size",
        type=int,
        required=False,
        default=10000,
        help="Window size for sliding windows.",
    )
    parser_count.add_argument(
        "-s",
        "--window-interval",
        type=int,
        required=False,
        default=10000,
        help="Window interval for sliding windows.",
    )
    parser_count.add_argument(
        "-k",
        "--kmer-size",
        type=int,
        required=False,
        default=3,
        help="Kmer size to count.",
    )
    parser_count.add_argument(
        "-t",
        "--threads",
        type=int,
        required=False,
        default=1,
        help="Number of threads for parallel processing.",
    )

    parser_count.set_defaults(func=run_count)

    # Subparser for plot
    parser_plot = subparsers.add_parser(
        "plot",
        help="Plot k-mer distributions and sequence signatures.",
        parents=[parent_parser],
    )
    parser_plot.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Path to input BEG file or database.",
    )
    parser_plot.add_argument(
        "-f",
        "--input-format",
        type=str,
        required=False,
        default="beg",
        choices=["beg", "db"],
        help="Input file format. Options: tab-separated values (beg); sqlite database (db).",
    )
    parser_plot.add_argument(
        "-o",
        "--output-prefix",
        type=str,
        required=False,
        default=None,
        help="Prefix for output images/figures.",
    )
    parser_plot.add_argument(
        "-p",
        "--output-format",
        type=str,
        required=False,
        default="png",
        choices=["png", "html"],
    )
    parser_plot.set_defaults(func=run_plot)

    return parser


def main():
    """Main CLI entry point."""
    parser = build_parser()
    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(0)

    configure_logging(getattr(args, "verbose", False))

    args.func(args)


if __name__ == "__main__":
    main()
