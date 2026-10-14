"""Command handler for sequence plotting."""

import logging

import src.beglib as bl  # noqa: F401

logger = logging.getLogger(__name__)


def run_plot(args) -> None:
    """Execute the plot command.

    Args:
        args: Parsed command line arguments containing input and output file paths.
    """
    logger.info("Executing 'plot' subparser command...")
