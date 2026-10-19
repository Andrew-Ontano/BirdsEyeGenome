# Birds Eye Genome (BEG)
*"Genomic sequence signature at a distance."*

**Birds Eye Genome** (**BEG**) is the successor to the [Spectra](https://github.com/Andrew-Ontano/Spectra) toolkit for visualizing genome-scale sequences and inferring regional shifts in coherent genome patterns. By supplying a FASTA file of one or more genomic sequences, **BEG** outputs both static figures for publication and interactive Web application modules for exploring genome-wide k-mer signatures and segmenting regional domain identities.

---

## Installation

### Prerequisites & Conda Environment Setup

BEG requires Python 3.12 or higher. We recommend using Conda or Mamba to create an isolated environment with all necessary scientific dependencies:

```bash
# Create a new Conda environment with Python 3.12
conda create -n beg_env python=3.12 -y

# Activate the environment
conda activate beg_env

# Clone repository and install BEG in editable mode
git clone https://github.com/Andrew-Ontano/Birds-Eye-Genome.git
cd Birds-Eye-Genome
pip install -e .
```

Dependencies installed via `pyproject.toml` include `biopython`, `pandas`, `numpy`, `matplotlib`, `plotly`, `dash`, `scikit-learn`, `scipy`, and `hmmlearn`.

---

## Workflow Synopsis

BEG provides a modular CLI interface via `beg.py` and dedicated utilities in `src/utils/`:

1. **`beg count`**: Sliding window k-mer frequency counting across FASTA sequence contigs. Computes k-mer counts alongside zlib compression complexity per window and outputs tab-separated `.beg` text files or binary SQLite `.db` databases.
2. **`beg domains`**: Dimensionality reduction (PCA) and sequence segmentation using Hidden Markov Models (GaussianHMM). Uncovers hidden domain states across expected cluster counts ($N$-states) and computes Shannon entropy and zlib complexity per window.
3. **`beg plot`**: Generates high-resolution stacked sequence signature figures (PNG, PDF, SVG) displaying k-mer composition patterns.
4. **`beg server`**: Launches an interactive Dash web application (`http://127.0.0.1:8050`) for interactive genome signature inspection.
5. **`src/utils/domain_visualize.py`**: Visualizes aligned stacked k-mer signatures, Shannon entropy line graphs, zlib complexity line graphs, and horizontal domain state prediction blocks.
6. **`src/utils/hierarchical_composition.py`**: Integrates multi-window domain datasets across different window sizes and state counts into a consolidated hierarchical genome composition table.

---

## Example Pipeline Block

A complete end-to-end multi-window analytical pipeline for a genome FASTA file (`genome.fasta`):

```bash
# 1. Count k-mers across small (10kb) and large (50kb) window sizes
beg count -i genome.fasta -w 10000 -s 10000 -o genome_10k.beg
beg count -i genome.fasta -w 50000 -s 50000 -o genome_50k.beg

# 2. Perform regional domain segmentation across state ranges (e.g. 2 to 4 states)
beg domains -i genome_10k.beg -n "2:4" -o genome_10k_domains.beg
beg domains -i genome_50k.beg -n "2:4" -o genome_50k_domains.beg

# 3. Generate static signature plots
beg plot -i genome_10k.beg -o genome_10k_signature --format png

# 4. Construct multi-window hierarchical composition
python src/utils/hierarchical_composition.py \
  -i genome_10k_domains.beg genome_50k_domains.beg \
  -o genome_hierarchical_composition.tsv

# 5. Render aligned domain visualization figure
python src/utils/domain_visualize.py \
  -d genome_10k_domains.beg \
  -c genome_10k.beg \
  -o figures/
```

---

## Example Calls to Utilities

### `domain_visualize.py`

Generate aligned sequence signature and HMM domain state figures:

```bash
python src/utils/domain_visualize.py \
  -d genome_10k_domains.beg \
  -c genome_10k.beg \
  -o output_figures/ \
  --format png \
  --dpi 300
```

### `hierarchical_composition.py`

Consolidate multi-window domain outputs into hierarchical labels and sequential levels:

```bash
python src/utils/hierarchical_composition.py \
  -i domain_50k.beg domain_25k.beg domain_10k.beg \
  -o hierarchical_output.tsv
```

---

## Information and Statistical Analyses Used

BEG integrates several statistical and algorithmic methods to quantify sequence signatures and partition regional genomic domain architectures:

1. **Sliding Window K-mer Composition**:
   - Computes dense distributions of short k-mers (e.g. trinucleotides, $k=3$, 64 possible combinations or 32 canonical reverse-complement pairs) across sliding genomic windows.

2. **zlib Compression Model Complexity**:
   - Measures sequence complexity using lossless zlib compression ratio ($C = \frac{\text{bytes}(\text{zlib}(\text{seq}))}{\text{len}(\text{seq})}$).
   - *Citation*: Compression ratios provide robust estimates of algorithmic sequence complexity and entropy rate (e.g., Li & Vitányi, *An Introduction to Kolmogorov Complexity and Its Applications*, Springer).

3. **Shannon Information Entropy**:
   - Calculates relative k-mer proportion entropy across window features using the natural logarithm ($H = -\sum_{i} p_i \ln(p_i)$).
   - *Citation*: Shannon, C. E. (1948). A mathematical theory of communication. *Bell System Technical Journal*, 27(3), 379-423.

4. **1D Gaussian Smoothing Filter**:
   - Applies Gaussian smoothing ($\sigma=1.0$) across adjacent spatial windows per sequence contig prior to PCA to reduce single-window noise while preserving boundaries.

5. **Principal Component Analysis (PCA)**:
   - Reduces high-dimensional normalized k-mer frequency matrices down to principal components capturing target cumulative variance ($\ge 90\%$).
   - *Citation*: Pearson, K. (1901). On lines and planes of closest fit to systems of points in space. *Philosophical Magazine*, 2(11), 559-572; Pedregosa et al. (2011). Scikit-learn: Machine Learning in Python. *JMLR*, 12, 2825-2830.

6. **Hidden Markov Model (HMM) Regional Segmentation**:
   - Models transitions between hidden domain states across PCA component representations using Gaussian HMMs with full covariance structure (`hmmlearn.hmm.GaussianHMM`).
   - *Citation*: Rabiner, L. R. (1989). A tutorial on Hidden Markov Models and selected applications in speech recognition. *Proceedings of the IEEE*, 77(2), 257-286.
