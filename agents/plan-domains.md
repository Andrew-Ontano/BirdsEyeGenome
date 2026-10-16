# Status: Completed
# Objective:
Develop a generalized Python pipeline to identify and classify genomic regional domains based on windowed k-mer composition. The pipeline must handle k-mers of any length, perform dimensionality reduction, support both full and canonical k-mer sets, and output labeled regional domains with their underlying k-mer frequencies.

## 1) Core Data Structures & Inputs
* **Input Data**: A Pandas DataFrame containing sequential genomic windows. Columns include metadata (`seq_name`, `window_start`, `window_stop`) and k-mer counts.
* **Target Output**: A detailed BEG/db where each window is assigned a `Domain_ID` and `Domain_Label`, alongside a summary table exporting the mean k-mer frequencies for each distinct domain type.
## 2) Pipeline Modules
### A. Preprocessing & K-mer canonicalization
* **Requirement**: Dynamically handle any k-mer size ($k$).
* **Canonical Collapsing**: Implement a function to generate canonical k-mers (the lexicographically smaller sequence between a k-mer and its reverse complement).
* **Matrix Transformation**:
* * _Full Set Pathway_: Convert the raw k-mer count matrix to row-wise proportions (frequencies).
* * _Canonical Set Pathway_: Sum the counts of each k-mer and its reverse complement, drop the redundant columns, and then convert to row-wise proportions.
### B. Dimensionality Reduction (PCA)
* **Requirement**: Prevent the segmentation algorithm from being overwhelmed by the $4^k$ feature space.
* **Scaling**: Apply `sklearn.preprocessing.StandardScaler` to the frequency matrices.
* **PCA Fitting**: Apply `sklearn.decomposition.PCA`.
* **Parameter**: Allow `n_components` to be dynamic, defaulting to capturing ~90% of the variance (usually 5-10 components for $k \le 6$).
### C. Domain Segmentation (Gaussian HMM)
* **Requirement**: Model the sequence of windows as hidden states emitted by underlying regional identities. This avoids rigid breakpoint penalties and handles variable chromosome lengths automatically.
* **Implementation**: Use `hmmlearn.hmm.GaussianHMM`.
* **Parameters**:
* * `n_components`: Number of expected domain identities (e.g., 3 to 7).
* * `covariance_type="full"`.
* **Execution**: Fit the HMM on the PCA-reduced space. Use `model.predict()` to assign a hidden state (`Domain ID`) to each contiguous genomic window.
### D. Aggregation & Labeling
* **Requirement**: Map the abstract HMM states back to biological meaning based on their k-mer profiles.
* **Aggregation**: Group the original DataFrame by the assigned Domain_ID. Calculate the mean frequency of every k-mer within each domain.
* **Labeling Logic**: Assign a `Domain_Label` based on the dominant compositional features of that state (e.g., "AT-rich", "GC-rich", "High-complexity"). Note to Jules: Implement a customizable scoring function here so the user can define labeling rules based on their specific k-mers of interest after the fact.
## 3) Execution Flow (Draft Python Logic for Jules)
```python
# Pseudo-code architecture for Jules to expand
def identify_regional_domains(
    df: pd.DataFrame, 
    kmer_cols: list[str], 
    use_canonical: bool = False, 
    n_states: int = 4, 
    variance_target: float = 0.90
) -> tuple[pd.DataFrame, pd.DataFrame]:
    
    # 1. Prepare frequency matrix
    if use_canonical:
        # Merge reverse complements
        df_counts = collapse_to_canonical(df, kmer_cols)
        canonical_cols = get_canonical_list(kmer_cols)
        freq_matrix = normalize_counts(df_counts[canonical_cols])
    else:
        freq_matrix = normalize_counts(df[kmer_cols])
        
    # 2. Dimensionality Reduction
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(freq_matrix)
    
    pca = PCA(n_components=variance_target)
    pca_components = pca.fit_transform(scaled_data)
    
    # 3. HMM Segmentation
    hmm_model = GaussianHMM(n_components=n_states, covariance_type="full", random_state=42)
    hmm_model.fit(pca_components)
    domain_assignments = hmm_model.predict(pca_components)
    
    # 4. Attach assignments back to original data
    df_results = df.copy()
    df_results["Domain_ID"] = domain_assignments
    
    # 5. Aggregate domain profiles
    domain_profiles = aggregate_domain_frequencies(df_results, kmer_cols, "Domain_ID")
    
    # 6. Apply Labeling (Extensible)
    df_results = assign_domain_labels(df_results, domain_profiles)
    
    return df_results, domain_profiles
```
## 4) Directives for Jules
* **Memory Efficiency**: Use vectorized Pandas and NumPy operations for the canonical collapsing to avoid fragmentation or slow iteration.
* **Missing Data**: Ensure the HMM cleanly handles or ignores empty windows (all-zero counts) without throwing NaN errors during scaling. N-gaps should not invalidate block identity.
* **Modularity**: Keep the canonical conversion, PCA, and HMM steps as separate, testable functions so they can be benchmarked independently.