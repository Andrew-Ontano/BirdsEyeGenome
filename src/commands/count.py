"""Command handler for k-mer counting."""

import logging
import sys

import src.beglib as bl

logger = logging.getLogger(__name__)


def run_count(args) -> None:
    """Execute the count command.

    Args:
        args: Parsed command line arguments containing input and output file paths.
    """
    logger.info("Executing 'count' subparser command...")
    try:
        records = bl.read_genome(args.input, args.input_format)
        record_count = 0
        for record in records:
            record_count += 1
            logger.info("Loaded accession: %s (length: %d bp)", record.id, len(record.seq))
        logger.info("Successfully loaded %d sequence record(s) from %s", record_count, args.input)
    except FileNotFoundError as err:
        logger.error("%s", err)
        sys.exit(1)
    except Exception as err:
        logger.error("Failed to parse genome file '%s': %s", args.input, err)
        sys.exit(1)
