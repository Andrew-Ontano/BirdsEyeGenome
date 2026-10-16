"""Command handler for regional domain segmentation and identification."""

import logging
import sqlite3
import sys
from pathlib import Path

import src.beglib as bl

logger = logging.getLogger(__name__)


def run_domains(args) -> None:
    """Execute the domains command to segment genomic windows into regional domains.

    Args:
        args: Parsed command line arguments containing input, output, and HMM settings.
    """
    logger.info("Executing 'domains' command...")
    try:
        input_path = Path(args.input)
        if not input_path.exists():
            logger.error("Input file not found: %s", args.input)
            sys.exit(1)

        input_format = getattr(args, "input_format", "auto")
        df = bl.load_counts_data(input_path, file_format=input_format)

        if df.empty:
            logger.error("Input dataset is empty: %s", args.input)
            sys.exit(1)

        meta_cols = ["library_name", "seq_name", "window_start", "window_stop", "window_size"]
        kmer_cols = [c for c in df.columns if c not in meta_cols]

        if not kmer_cols:
            logger.error("No k-mer count columns found in dataset.")
            sys.exit(1)

        use_canonical = getattr(args, "canonical", False)
        n_states = getattr(args, "n_states", 4)
        variance_target = getattr(args, "variance_target", 0.90)

        logger.info(
            "Running domain identification: n_states=%d, variance_target=%.2f, canonical=%s",
            n_states,
            variance_target,
            use_canonical,
        )

        df_results, domain_profiles = bl.identify_regional_domains(
            df=df,
            kmer_cols=kmer_cols,
            use_canonical=use_canonical,
            n_states=n_states,
            variance_target=variance_target,
        )

        output_format = getattr(args, "output_format", "auto")
        if output_format == "auto":
            output_format = "db" if input_path.suffix.lower() == ".db" else "beg"

        db_output_flag = getattr(args, "db_output", False)

        # Main window output path
        output_file = args.output
        if not output_file:
            ext = ".db" if output_format == "db" else ".beg"
            output_file = str(input_path.with_name(f"{input_path.stem}_domains{ext}"))

        summary_file = args.summary_output
        if not summary_file:
            summary_file = str(input_path.with_name(f"{input_path.stem}_domain_profiles.tsv"))

        # Save main window output (TSV .beg or SQLite .db)
        if output_format == "db":
            logger.info("Writing domain window results and summary profiles to SQLite DB: %s", output_file)
            run_id = f"run_n{n_states}_v{variance_target}_c{int(use_canonical)}"
            conn = sqlite3.connect(output_file)
            params = {
                "n_states": n_states,
                "variance_target": variance_target,
                "canonical": use_canonical,
                "input_file": str(input_path.resolve()),
            }
            bl.write_domains_db(conn, run_id, df_results, domain_profiles, params=params)
            conn.close()
        else:
            logger.info("Writing domain window results to BEG TSV file: %s", output_file)
            df_results.to_csv(output_file, sep="\t", index=False)

        # Save domain profiles summary TSV
        logger.info("Writing domain profiles summary to TSV file: %s", summary_file)
        domain_profiles.to_csv(summary_file, sep="\t", index=False)

        # If --db-output flag supplied when input is not DB or distinct DB requested
        if db_output_flag and output_format != "db":
            db_file = str(input_path.with_name(f"{input_path.stem}_domains.db"))
            logger.info("Writing domain run to supplementary SQLite DB: %s", db_file)
            run_id = f"run_n{n_states}_v{variance_target}_c{int(use_canonical)}"
            conn = sqlite3.connect(db_file)
            params = {
                "n_states": n_states,
                "variance_target": variance_target,
                "canonical": use_canonical,
                "input_file": str(input_path.resolve()),
            }
            bl.write_domains_db(conn, run_id, df_results, domain_profiles, params=params)
            conn.close()

        logger.info("Successfully identified regional domains for %d windows.", len(df_results))

    except Exception as err:
        logger.error("Failed to run domain identification: %s", err)
        sys.exit(1)
