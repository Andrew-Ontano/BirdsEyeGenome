"""Birds Eye Genome core library module."""

from __future__ import annotations

import itertools
import json
import sqlite3
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import Any

import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def read_genome(
    input_file: str,
    input_format: str = "fasta",
) -> Iterator[SeqRecord]:
    """Check that the input file exists and load sequence records via Biopython SeqIO.

    Args:
        input_file (str): Path to the input genomic sequence file.
        input_format (str, optional): Input sequence file format recognized by SeqIO.
            Defaults to "fasta".

    Returns:
        Iterator[SeqRecord]: An iterator yielding SeqRecord objects from the input file.

    Raises:
        FileNotFoundError: If the specified file does not exist or is not a regular file.
    """
    file_path = Path(input_file)
    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(f"Input file does not exist or is not a file: {input_file}")

    return SeqIO.parse(str(file_path), input_format)


def set_mers(kmer_size: int = 3) -> list[str]:
    """Build all possible k-mers constructed from A, C, G, and T.

    Args:
        kmer_size (int, optional): Length of the k-mers to generate. Defaults to 3.

    Returns:
        list[str]: Sorted list of all k-mer strings constructed from A, C, G, T.
    """
    nucleotides = ["A", "C", "G", "T"]
    kmers = ["".join(p) for p in itertools.product(nucleotides, repeat=kmer_size)]
    return sorted(kmers)


def window_task_generator(
    records: Iterable[SeqRecord],
    window_size: int,
    window_interval: int,
    kmer_list: list[str],
) -> Iterator[dict[str, Any]]:
    """Define sliding windows for each genomic sequence and yield window tasks.

    Args:
        records (Iterable[SeqRecord]): Iterable yielding SeqRecord objects.
        window_size (int): The window size in base pairs.
        window_interval (int): The interval step size between sliding window start positions.
        kmer_list (list[str]): The list of k-mers to be counted.

    Yields:
        dict[str, Any]: A dictionary containing window task parameters:
            - "seq_name": Sequence accession/ID.
            - "start": Start position (0-based inclusive).
            - "end": End position (0-based exclusive).
            - "window_size": Actual window length in base pairs.
            - "sequence": Subsequence string for the window.
            - "kmers": List of k-mer strings to count.
    """
    for record in records:
        seq_str = str(record.seq).upper()
        seq_len = len(seq_str)
        if seq_len == 0:
            continue

        for start in range(0, seq_len, window_interval):
            end = min(start + window_size, seq_len)
            subseq = seq_str[start:end]
            actual_window_size = len(subseq)
            yield {
                "seq_name": str(record.id),
                "start": start,
                "end": end,
                "window_size": actual_window_size,
                "sequence": subseq,
                "kmers": kmer_list,
            }


def count_mers(task: dict[str, Any]) -> dict[str, Any]:
    """Count occurrences of each k-mer within a window sequence task.

    Args:
        task (dict[str, Any]): Window task dictionary containing "seq_name", "start",
            "end", "window_size", "sequence", and "kmers".

    Returns:
        dict[str, Any]: Dictionary containing task metadata and a "counts" mapping of
            each k-mer to its count in the window.
    """
    subseq = task["sequence"].upper()
    kmers = task["kmers"]
    kmer_len = len(kmers[0]) if kmers else 0

    counts = {kmer: 0 for kmer in kmers}
    seq_len = len(subseq)

    if kmer_len > 0 and seq_len >= kmer_len:
        for i in range(seq_len - kmer_len + 1):
            kmer = subseq[i : i + kmer_len]
            if kmer in counts:
                counts[kmer] += 1

    return {
        "seq_name": task["seq_name"],
        "start": task["start"],
        "end": task["end"],
        "window_size": task["window_size"],
        "counts": counts,
    }


def init_db(db_path: str, kmers: list[str] | None = None, kmer_size: int = 3) -> sqlite3.Connection:
    db_file = Path(db_path)
    if db_file.exists():
        db_file.unlink()

    conn = sqlite3.connect(str(db_file))
    cursor = conn.cursor()

    cursor.execute("""
    CREATE TABLE IF NOT EXISTS metadata (
        key TEXT PRIMARY KEY,
        value TEXT
    )
    """)

    # Added the 'counts' BLOB column
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS windows (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        library_name TEXT,
        seq_name TEXT,
        window_start INTEGER,
        window_stop INTEGER,
        window_size INTEGER,
        counts BLOB 
    )
    """)

    cursor.execute("CREATE INDEX IF NOT EXISTS idx_windows_seq_name ON windows(seq_name)")

    if kmers and len(kmers) > 0:
        inferred_kmer_size = len(kmers[0])
        # Save the exact order of kmers so we can decode the binary array later
        cursor.execute("INSERT INTO metadata (key, value) VALUES (?, ?)", ("kmers_list", json.dumps(kmers)))
    else:
        inferred_kmer_size = kmer_size

    cursor.execute("INSERT OR REPLACE INTO metadata (key, value) VALUES (?, ?)", ("kmer_size", str(inferred_kmer_size)))
    conn.commit()
    return conn


def write_counts_db(
        conn: sqlite3.Connection,
        library_name: str,
        results: list[dict[str, Any]],
        kmers: list[str],
        kmer_size: int | None = None,
) -> None:
    cursor = conn.cursor()

    if kmer_size is None and kmers and len(kmers) > 0:
        kmer_size = len(kmers[0])
    elif kmer_size is None:
        kmer_size = 3

    cursor.execute("INSERT OR REPLACE INTO metadata (key, value) VALUES (?, ?)", ("kmer_size", str(kmer_size)))
    cursor.execute("INSERT OR REPLACE INTO metadata (key, value) VALUES (?, ?)", ("library_name", library_name))

    # Create a quick lookup map to place counts in the correct array index
    kmer_to_idx = {kmer: i for i, kmer in enumerate(kmers)}
    num_kmers = len(kmers)

    window_rows = []

    for res in results:
        # Initialize a dense array of zeros
        # uint32 allows counts up to 4.2 billion; use uint16 if max count < 65,535 for half the storage
        dense_counts = np.zeros(num_kmers, dtype=np.uint32)

        counts_dict = res.get("counts", {})
        for kmer, count in counts_dict.items():
            if count > 0 and kmer in kmer_to_idx:
                dense_counts[kmer_to_idx[kmer]] = count

        window_rows.append((
            library_name,
            res["seq_name"],
            res["start"] + 1,
            res["end"],
            res["window_size"],
            dense_counts.tobytes()  # Convert the numpy array directly to binary
        ))

    if window_rows:
        # One single executemany handles everything
        cursor.executemany(
            "INSERT INTO windows (library_name, seq_name, window_start, window_stop, window_size, counts) VALUES (?, ?, ?, ?, ?, ?)",
            window_rows,
        )

    conn.commit()


def write_counts_beg(
    output_file: str,
    input_file: str,
    library_name: str,
    kmers: list[str],
    results: list[dict[str, Any]],
) -> None:
    """Write window count results to a tab-separated BEG file with comment metadata.

    Args:
        output_file (str): Path to destination BEG file.
        input_file (str): Total path to the input sequence file.
        library_name (str): Shortened library name.
        kmers (list[str]): List of k-mers.
        results (list[dict[str, Any]]): List of window count result dictionaries.
    """
    abs_input_path = Path(input_file).resolve()
    with open(output_file, "w", encoding="utf-8") as f:
        f.write(f"# Input file: {abs_input_path}\n")
        f.write("# Generated by Birds Eye Genome (beg count)\n")

        header = ["library_name", "seq_name", "window_start", "window_stop", "window_size"] + kmers
        f.write("\t".join(header) + "\n")

        for res in results:
            counts = res.get("counts", {})
            row = [
                library_name,
                res["seq_name"],
                str(res["start"] + 1),
                str(res["end"]),
                str(res["window_size"]),
            ] + [str(counts.get(kmer, 0)) for kmer in kmers]
            f.write("\t".join(row) + "\n")


def kmer_sort_key(kmer: str) -> tuple[int, tuple[int, ...]]:
    """Sorting key for k-mers: AT-rich first (AAA > TTT > CCC > GGG).

    Args:
        kmer (str): K-mer sequence string.

    Returns:
        tuple[int, tuple[int, ...]]: Sorting tuple key.
    """
    gc_count = sum(1 for char in kmer if char in ("G", "C"))
    base_order = {"A": 0, "T": 1, "C": 2, "G": 3}
    base_ranks = tuple(base_order.get(char, 9) for char in kmer)
    return (gc_count, base_ranks)


def generate_kmer_colors(kmers: list[str]) -> dict[str, str]:
    """Generate color scheme for k-mers based on nucleotide composition.

    A-rich k-mers are red, T-rich are orange, C-rich are blue, and G-rich are greenish/teal.

    Args:
        kmers (list[str]): List of k-mer strings.

    Returns:
        dict[str, str]: Map of k-mer string to hex color code.
    """
    if not kmers:
        return {}

    k = len(kmers[0]) if kmers else 1

    base_colors = {
        "A": np.array([0.90, 0.15, 0.15]),
        "T": np.array([0.95, 0.50, 0.15]),
        "C": np.array([0.15, 0.40, 0.90]),
        "G": np.array([0.15, 0.75, 0.55]),
    }

    comp_groups: dict[tuple[int, int, int, int], list[str]] = {}
    for kmer in kmers:
        cA = kmer.count("A")
        cT = kmer.count("T")
        cC = kmer.count("C")
        cG = kmer.count("G")
        comp = (cA, cT, cC, cG)
        comp_groups.setdefault(comp, []).append(kmer)

    kmer_colors = {}
    for comp, kmer_list in comp_groups.items():
        cA, cT, cC, cG = comp
        fA, fT, fC, fG = cA / k, cT / k, cC / k, cG / k

        rgb_blend = fA * base_colors["A"] + fT * base_colors["T"] + fC * base_colors["C"] + fG * base_colors["G"]
        base_hsv = mcolors.rgb_to_hsv(rgb_blend)

        kmer_list.sort(key=kmer_sort_key)
        n = len(kmer_list)
        if n > 1:
            lum_steps = np.linspace(max(0.3, base_hsv[2] - 0.2), min(0.95, base_hsv[2] + 0.2), n)
        else:
            lum_steps = [base_hsv[2]]

        for i, kmer in enumerate(kmer_list):
            new_hsv = (base_hsv[0], base_hsv[1], lum_steps[i])
            new_rgb = mcolors.hsv_to_rgb(new_hsv)
            kmer_colors[kmer] = mcolors.to_hex(new_rgb)

    return kmer_colors


def get_reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence string.

    Args:
        seq (str): Nucleotide sequence string.

    Returns:
        str: Reverse complemented nucleotide sequence.
    """
    trans = str.maketrans("ATCGatcg", "TAGCtagc")
    return seq.translate(trans)[::-1]


def get_canonical_kmer(kmer: str) -> str:
    """Return the lexicographically smaller k-mer between kmer and its reverse complement.

    Args:
        kmer (str): K-mer string.

    Returns:
        str: Lexicographically smaller canonical k-mer.
    """
    rc = get_reverse_complement(kmer)
    return min(kmer, rc)


def load_counts_data(file_path: str, file_format: str = "auto") -> pd.DataFrame:
    """Load k-mer count data from a .beg TSV file or .db SQLite database into a DataFrame.

    Args:
        file_path (str): Path to the .beg or .db file.
        file_format (str): "beg", "db", or "auto" to infer from file extension.

    Returns:
        pd.DataFrame: DataFrame containing window counts.
    """
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"Input file does not exist: {file_path}")

    fmt = file_format.lower()
    if fmt == "auto":
        fmt = "db" if path.suffix.lower() == ".db" else "beg"

    if fmt == "db":
        conn = sqlite3.connect(str(path))
        cursor = conn.cursor()

        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='windows'")
        has_windows_table = cursor.fetchone() is not None

        if has_windows_table:
            # Check if using the new BLOB schema by looking for the 'counts' column
            cursor.execute("PRAGMA table_info(windows)")
            columns = [col[1] for col in cursor.fetchall()]
            is_dense_schema = "counts" in columns

            if is_dense_schema:
                # --- NEW DENSE BLOB SCHEMA ---
                try:
                    cursor.execute("SELECT value FROM metadata WHERE key = 'kmers_list'")
                    row = cursor.fetchone()
                    if row and row[0]:
                        kmers = json.loads(row[0])
                    else:
                        raise ValueError("Dense schema requires 'kmers_list' in metadata.")
                except sqlite3.OperationalError:
                    raise ValueError("Metadata table missing or corrupted.")

                # Fetch window metadata and the counts BLOB
                cursor.execute(
                    "SELECT library_name, seq_name, window_start, window_stop, window_size, counts "
                    "FROM windows ORDER BY id ASC"
                )
                rows = cursor.fetchall()
                conn.close()

                if not rows:
                    return pd.DataFrame()

                # Pre-allocate lists for metadata and arrays for counts
                meta_cols = {"library_name": [], "seq_name": [], "window_start": [], "window_stop": [],
                             "window_size": []}
                counts_list = []

                for row in rows:
                    meta_cols["library_name"].append(row[0])
                    meta_cols["seq_name"].append(row[1])
                    meta_cols["window_start"].append(row[2])
                    meta_cols["window_stop"].append(row[3])
                    meta_cols["window_size"].append(row[4])

                    if row[5]:
                        counts_list.append(np.frombuffer(row[5], dtype=np.uint32))
                    else:
                        counts_list.append(np.zeros(len(kmers), dtype=np.uint32))

                # Build the dataframe from the vectorized array to prevent fragmentation
                df_meta = pd.DataFrame(meta_cols)
                df_counts = pd.DataFrame(np.vstack(counts_list), columns=kmers)

                # Combine metadata and counts block
                df = pd.concat([df_meta, df_counts], axis=1)

            else:
                # --- OLD SPARSE EAV SCHEMA (Backwards compatibility) ---
                kmer_size = 3
                try:
                    cursor.execute("SELECT value FROM metadata WHERE key = 'kmer_size'")
                    row = cursor.fetchone()
                    if row and row[0]:
                        kmer_size = int(row[0])
                except sqlite3.OperationalError:
                    pass

                windows_df = pd.read_sql_query(
                    "SELECT id, library_name, seq_name, window_start, window_stop, window_size FROM windows ORDER BY id ASC",
                    conn,
                )
                counts_df = pd.read_sql_query("SELECT window_id, kmer, count FROM kmer_counts", conn)
                conn.close()

                if windows_df.empty:
                    return pd.DataFrame()

                if not counts_df.empty and "kmer" in counts_df.columns:
                    non_empty_kmers = counts_df["kmer"].dropna()
                    if not non_empty_kmers.empty:
                        kmer_size = len(str(non_empty_kmers.iloc[0]))

                # Note: Assumes `set_mers` is imported/available in your file scope
                all_kmers = set_mers(kmer_size)

                if not counts_df.empty:
                    pivoted = counts_df.pivot(index="window_id", columns="kmer", values="count").fillna(0).astype(int)
                else:
                    pivoted = pd.DataFrame(index=windows_df["id"])

                pivoted = pivoted.reindex(columns=all_kmers, fill_value=0)
                df = windows_df.merge(pivoted, left_on="id", right_index=True, how="left")
                df = df.fillna(0)
                if "id" in df.columns:
                    df = df.drop(columns=["id"])
        else:
            # Absolute fallback
            df = pd.read_sql_query("SELECT * FROM kmer_counts ORDER BY id ASC", conn)
            conn.close()
            if "id" in df.columns:
                df = df.drop(columns=["id"])
    else:
        # TSV Path
        df = pd.read_csv(str(path), sep="\t", comment="#")

    return df


def process_window_proportions(
    df: pd.DataFrame,
    canonical: bool = False,
    show_n_gaps: bool = False,
) -> tuple[pd.DataFrame, list[str], dict[str, str]]:
    """Process raw k-mer window counts into stacked proportion values for plotting.

    Args:
        df (pd.DataFrame): Input DataFrame containing window count data.
        canonical (bool): If True, collapse reverse-complement k-mers into canonical k-mers.
        show_n_gaps (bool): If True, calculate and include N-gaps proportion layer.

    Returns:
        tuple[pd.DataFrame, list[str], dict[str, str]]:
            - Processed DataFrame with proportion values.
            - List of feature column names (in order of stacked display).
            - Map of feature name to hex color string.
    """
    meta_cols = ["library_name", "seq_name", "window_start", "window_stop", "window_size"]
    kmer_cols = [c for c in df.columns if c not in meta_cols]

    if not kmer_cols:
        raise ValueError("No k-mer count columns found in the dataset.")

    processed_df = df[meta_cols].copy()
    counts_df = df[kmer_cols].copy()

    kmer_len = len(kmer_cols[0])

    if canonical:
        canonical_map: dict[str, list[str]] = {}
        for kmer in kmer_cols:
            canon = get_canonical_kmer(kmer)
            canonical_map.setdefault(canon, []).append(kmer)

        canon_counts = pd.DataFrame(index=df.index)
        for canon in sorted(canonical_map.keys(), key=kmer_sort_key):
            cols_to_sum = canonical_map[canon]
            canon_counts[canon] = counts_df[cols_to_sum].sum(axis=1)

        counts_df = canon_counts
        feature_kmers = list(counts_df.columns)
    else:
        feature_kmers = sorted(kmer_cols, key=kmer_sort_key)
        counts_df = counts_df[feature_kmers]

    raw_kmer_sum = df[kmer_cols].sum(axis=1)

    feature_cols = []
    color_map = {}

    if show_n_gaps:
        ns_count = (df["window_size"] - raw_kmer_sum - kmer_len + 1).clip(lower=0)
        counts_df["Ns"] = ns_count
        feature_cols.append("Ns")
        color_map["Ns"] = "#000000"

    feature_cols.extend(feature_kmers)

    kmer_color_mapping = generate_kmer_colors(feature_kmers)
    color_map.update(kmer_color_mapping)

    total_counts = counts_df[feature_cols].sum(axis=1)
    norm_denom = total_counts.replace(0, 1)

    normalized_features = counts_df[feature_cols].div(norm_denom, axis=0)
    processed_df = pd.concat([processed_df, normalized_features], axis=1)

    return processed_df, feature_cols, color_map


def collapse_to_canonical(df: pd.DataFrame, kmer_cols: list[str]) -> pd.DataFrame:
    """Collapse k-mer count columns into canonical k-mers (summing reverse complements).

    Args:
        df (pd.DataFrame): DataFrame containing k-mer counts.
        kmer_cols (list[str]): List of k-mer column names to collapse.

    Returns:
        pd.DataFrame: DataFrame with redundant reverse complement k-mers merged into
            lexicographically smaller canonical k-mer columns.
    """
    canonical_map: dict[str, list[str]] = {}
    for kmer in kmer_cols:
        canon = get_canonical_kmer(kmer)
        canonical_map.setdefault(canon, []).append(kmer)

    df_res = pd.DataFrame(index=df.index)
    sorted_canons = sorted(canonical_map.keys(), key=kmer_sort_key)
    for canon in sorted_canons:
        cols_to_sum = canonical_map[canon]
        df_res[canon] = df[cols_to_sum].sum(axis=1)

    return df_res


def prepare_kmer_frequencies(
    df: pd.DataFrame,
    kmer_cols: list[str],
    use_canonical: bool = False,
) -> tuple[pd.DataFrame, list[str]]:
    """Convert raw window k-mer counts to row-wise relative frequencies (proportions).

    Args:
        df (pd.DataFrame): Input DataFrame containing k-mer count columns.
        kmer_cols (list[str]): Names of k-mer count columns in df.
        use_canonical (bool): If True, collapse counts into canonical k-mers before normalizing.

    Returns:
        tuple[pd.DataFrame, list[str]]:
            - DataFrame of normalized frequency values (0.0 to 1.0).
            - List of resulting feature column names.
    """
    if use_canonical:
        counts_df = collapse_to_canonical(df, kmer_cols)
        feature_cols = list(counts_df.columns)
    else:
        feature_cols = sorted(kmer_cols, key=kmer_sort_key)
        counts_df = df[feature_cols].copy()

    row_sums = counts_df.sum(axis=1)
    denom = row_sums.replace(0, 1.0)
    freq_df = counts_df.div(denom, axis=0)

    return freq_df, feature_cols


def parse_n_states(n_states_arg: str | int) -> list[int]:
    """Parse n_states argument into a list of integer state numbers.

    Supports singular integer values (e.g., "4" or 4) or colon-separated ranges
    (e.g., "3:10"), inclusive of limits.

    Args:
        n_states_arg (str | int): Singular integer or range string "A:B".

    Returns:
        list[int]: Sorted list of state numbers.
    """
    arg_str = str(n_states_arg).strip()
    if ":" in arg_str:
        parts = arg_str.split(":")
        if len(parts) != 2:
            raise ValueError(f"Invalid range format for n-states: '{n_states_arg}'. Expected 'start:stop'.")
        start, stop = int(parts[0]), int(parts[1])
        if start > stop:
            raise ValueError(f"Invalid range for n-states: '{n_states_arg}'. Start ({start}) > stop ({stop}).")
        return list(range(start, stop + 1))
    return [int(arg_str)]


def identify_regional_domains(
    df: pd.DataFrame,
    kmer_cols: list[str],
    use_canonical: bool = False,
    n_states: int | str = 4,
    variance_target: float = 0.90,
    random_state: int = 42,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Segment genomic windows into hidden Markov model domain states based on k-mer composition.

    Supports evaluating a single state count or iterating across a range of expected state values.

    Args:
        df (pd.DataFrame): Input genomic window DataFrame containing metadata and k-mer counts.
        kmer_cols (list[str]): List of k-mer count column names.
        use_canonical (bool): Collapse reverse complement k-mers if True. Defaults to False.
        n_states (int | str): Number of HMM states or range string "A:B" (e.g. "3:10").
            Defaults to 4.
        variance_target (float): Variance ratio target for PCA dimensionality reduction.
            Defaults to 0.90.
        random_state (int): Seed for HMM initialization reproducible results.

    Returns:
        tuple[pd.DataFrame, pd.DataFrame]:
            - Truncated window DataFrame with localization columns (`library_name`, `seq_name`,
              `window_start`, `window_stop`, `window_size`), `comment`, and `Expected_{N}` columns.
            - Combined domain profiles DataFrame containing mean k-mer frequencies for each domain.
    """
    from hmmlearn.hmm import GaussianHMM
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    if df.empty:
        raise ValueError("Input DataFrame is empty.")

    state_list = parse_n_states(n_states)

    meta_cols = ["library_name", "seq_name", "window_start", "window_stop", "window_size"]
    loc_cols = [c for c in meta_cols if c in df.columns]

    # Prepare truncated output DataFrame
    df_results = df[loc_cols].copy() if loc_cols else pd.DataFrame(index=df.index)

    # 1. Prepare k-mer relative frequencies
    freq_df, feature_cols = prepare_kmer_frequencies(df, kmer_cols, use_canonical=use_canonical)

    # Detect zero/empty windows (e.g. N-gaps)
    zero_mask = freq_df.sum(axis=1) == 0.0

    # Set comments for zero/N-gap windows
    df_results["comment"] = ["N-gap present" if zero_mask.iloc[i] else "" for i in range(len(df))]

    # 2. Handle missing/empty windows via forward/backward fill interpolation
    model_freq_df = freq_df.copy()
    if zero_mask.any():
        model_freq_df.loc[zero_mask, :] = np.nan
        model_freq_df = model_freq_df.ffill().bfill().fillna(0.0)

    # 3. Dimensionality Reduction (StandardScaler + PCA)
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(model_freq_df)

    max_components = min(scaled_data.shape[0], scaled_data.shape[1])
    if max_components < 1:
        raise ValueError("Not enough features or samples for PCA.")

    pca_target: float | int = variance_target
    if isinstance(variance_target, float) and (variance_target <= 0 or variance_target >= 1.0):
        pca_target = min(5, max_components)
    elif isinstance(variance_target, float) and max_components == 1:
        pca_target = 1

    pca = PCA(n_components=pca_target, random_state=random_state)
    pca_components = pca.fit_transform(scaled_data)

    lengths = df.groupby("seq_name", sort=False).size().tolist() if "seq_name" in df.columns else None

    profile_frames = []

    # Iterate over all requested state expectations
    for num_s in state_list:
        effective_states = min(num_s, pca_components.shape[0])
        hmm_model = GaussianHMM(
            n_components=effective_states,
            covariance_type="full",
            random_state=random_state,
            n_iter=100,
        )

        if lengths:
            hmm_model.fit(pca_components, lengths=lengths)
            domain_assignments = hmm_model.predict(pca_components, lengths=lengths)
        else:
            hmm_model.fit(pca_components)
            domain_assignments = hmm_model.predict(pca_components)

        col_name = f"Expected_{num_s}"
        df_results[col_name] = [f"Domain_{s}" for s in domain_assignments]

        # Aggregate summary profile for this expected state count
        temp_df = pd.DataFrame({"Domain_ID": df_results[col_name].values})
        prof = aggregate_domain_frequencies(temp_df, feature_cols, freq_df)
        prof.insert(0, "Expected_States", num_s)
        profile_frames.append(prof)

    combined_profiles = pd.concat(profile_frames, ignore_index=True)

    return df_results, combined_profiles


def aggregate_domain_frequencies(
    df_results: pd.DataFrame,
    feature_cols: list[str],
    freq_df: pd.DataFrame,
) -> pd.DataFrame:
    """Calculate mean k-mer frequencies for each distinct domain identity state.

    Args:
        df_results (pd.DataFrame): DataFrame containing `Domain_ID` column.
        feature_cols (list[str]): List of k-mer feature column names.
        freq_df (pd.DataFrame): DataFrame of k-mer frequencies aligned with df_results.

    Returns:
        pd.DataFrame: Summary table with Domain_ID as index/column and mean frequency per k-mer.
    """
    freq_with_domain = freq_df.copy()
    freq_with_domain["Domain_ID"] = df_results["Domain_ID"].values

    profiles = freq_with_domain.groupby("Domain_ID", as_index=False)[feature_cols].mean()
    return profiles


def write_domains_db(
    conn: sqlite3.Connection,
    run_id: str,
    df_results: pd.DataFrame,
    domain_profiles: pd.DataFrame,
    params: dict[str, Any] | None = None,
) -> None:
    """Write domain assignments and domain summary profiles to SQLite tables.

    Args:
        conn (sqlite3.Connection): Active SQLite connection.
        run_id (str): Unique identifier string for this run parameterization.
        df_results (pd.DataFrame): Window assignments containing seq_name, window_start,
            window_stop, comment, and Expected_{N} columns.
        domain_profiles (pd.DataFrame): Summary table of mean k-mer frequencies per domain.
        params (dict[str, Any], optional): Dictionary of parameter metadata to log.
    """
    cursor = conn.cursor()

    cursor.execute("""
    CREATE TABLE IF NOT EXISTS domain_runs (
        run_id TEXT PRIMARY KEY,
        created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
        params TEXT
    )
    """)

    cursor.execute("""
    CREATE TABLE IF NOT EXISTS domains (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        run_id TEXT,
        seq_name TEXT,
        window_start INTEGER,
        window_stop INTEGER,
        assignments_json TEXT,
        comment TEXT,
        FOREIGN KEY(run_id) REFERENCES domain_runs(run_id)
    )
    """)

    cursor.execute("""
    CREATE TABLE IF NOT EXISTS domain_profiles (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        run_id TEXT,
        expected_states INTEGER,
        domain_id TEXT,
        profile_json TEXT,
        FOREIGN KEY(run_id) REFERENCES domain_runs(run_id)
    )
    """)

    cursor.execute("CREATE INDEX IF NOT EXISTS idx_domains_run_seq ON domains(run_id, seq_name)")

    # Insert run metadata
    params_json = json.dumps(params or {})
    cursor.execute(
        "INSERT OR REPLACE INTO domain_runs (run_id, params) VALUES (?, ?)",
        (run_id, params_json),
    )

    # Clear previous domain results for this run_id if re-running
    cursor.execute("DELETE FROM domains WHERE run_id = ?", (run_id,))
    cursor.execute("DELETE FROM domain_profiles WHERE run_id = ?", (run_id,))

    expected_cols = [c for c in df_results.columns if c.startswith("Expected_")]

    # Insert window domain assignments
    domain_rows = []
    for _, row in df_results.iterrows():
        assign_dict = {col: str(row[col]) for col in expected_cols}
        domain_rows.append((
            run_id,
            str(row.get("seq_name", "")),
            int(row.get("window_start", 0)),
            int(row.get("window_stop", 0)),
            json.dumps(assign_dict),
            str(row.get("comment", "")),
        ))

    if domain_rows:
        cursor.executemany(
            "INSERT INTO domains (run_id, seq_name, window_start, window_stop, assignments_json, comment) "
            "VALUES (?, ?, ?, ?, ?, ?)",
            domain_rows,
        )

    # Insert domain profiles
    feature_cols = [c for c in domain_profiles.columns if c not in ("Expected_States", "Domain_ID")]
    profile_rows = []
    for _, row in domain_profiles.iterrows():
        exp_s = int(row.get("Expected_States", 0))
        domain_id = str(row["Domain_ID"])
        p_dict = {col: float(row[col]) for col in feature_cols}
        profile_rows.append((
            run_id,
            exp_s,
            domain_id,
            json.dumps(p_dict),
        ))

    if profile_rows:
        cursor.executemany(
            "INSERT INTO domain_profiles (run_id, expected_states, domain_id, profile_json) VALUES (?, ?, ?, ?)",
            profile_rows,
        )

    conn.commit()
