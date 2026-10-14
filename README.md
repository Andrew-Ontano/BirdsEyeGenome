# Birds Eye Genome
*"Genomic sequence signature at a distance."*

**Birds Eye Genome** (**BEG**) is the successor to the [Spectra](https://github.com/Andrew-Ontano/Spectra) toolkit for visualizing genome-scale sequences and inferring regional shifts in coherent genome patterns. By supplying a fasta of one or more genomic sequences, **BEG** outputs both static images ready for inclusion manuscript figures or interactive webpages for viewing large-scale DNA patterns.

**BEG** works via assessing the pattern of short k-mers (trinucleotides) and measuring their prevalence across each sequence via sliding windows. This results in a BEG file with the window positions and corresponding k-mer counts. BEG files are then used to generate plots or perform compositional analysis.