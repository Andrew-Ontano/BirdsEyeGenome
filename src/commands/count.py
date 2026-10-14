"""Command handler for k-mer counting."""

import logging
import multiprocessing
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
        kmers = bl.set_mers(args.kmer_size)
        logger.info("Generated %d k-mers for k-mer size %d", len(kmers), args.kmer_size)

        tasks = list(
            bl.window_task_generator(
                records,
                window_size=args.window_size,
                window_interval=args.window_interval,
                kmer_list=kmers,
            )
        )
        logger.info("Generated %d window task(s)", len(tasks))

        threads = getattr(args, "threads", 1)
        results = []

        if threads > 1 and len(tasks) > 0:
            logger.info("Processing window tasks using %d worker threads/processes", threads)
            with multiprocessing.Pool(processes=threads) as pool:
                results = pool.map(bl.count_mers, tasks)
        else:
            logger.info("Processing window tasks sequentially")
            results = [bl.count_mers(task) for task in tasks]

        logger.info("Successfully completed k-mer counting for %d window(s)", len(results))

    except FileNotFoundError as err:
        logger.error("%s", err)
        sys.exit(1)
    except Exception as err:
        logger.error("Failed to parse genome file '%s': %s", args.input, err)
        sys.exit(1)
