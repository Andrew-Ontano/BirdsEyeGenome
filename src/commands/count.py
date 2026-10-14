"""Command handler for k-mer counting."""

import argparse

import src.beglib as bl


def run_count(args: argparse.Namespace) -> None:
    """Execute the count command.

    Args:
        args: Parsed command line arguments containing input and output file paths.
    """
    print("Executing 'count' subparser command...")
    bl.count_kmers(input_path=args.input, output_path=args.output)
