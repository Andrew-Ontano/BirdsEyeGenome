"""Command handler for sequence plotting."""

import argparse

import src.beglib as bl


def run_plot(args: argparse.Namespace) -> None:
    """Execute the plot command.

    Args:
        args: Parsed command line arguments containing input and output file paths.
    """
    print("Executing 'plot' subparser command...")
    bl.plot(input_path=args.input, output_path=args.output)
