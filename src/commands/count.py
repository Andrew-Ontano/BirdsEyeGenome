"""Command handler for k-mer counting."""
import src.beglib as bl


def run_count(args):
    """Execute the count command.

    Args:
        args: Parsed command line arguments containing input and output file paths.
    """
    print("Executing 'count' subparser command...")