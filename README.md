Bioinformatics Assignment 4

Author: Avishka Kavindi Vithanage (s225407945)
Date: 2025-10-06

Overview

This repository contains my solutions for Bioinformatics Assignment 4, split into two separate R Markdown reports as requested:

Part 1 – Data analysis (RNA-seq counts + tree growth)

Part 2 – Sequence diversity (E. coli vs Cryobacterium sp. SO1)

Both reports are fully reproducible: they create needed folders, download inputs (where applicable), write tables to output/, and save figures to figures/.

Contents
bioinfo-assignment4/
├─ Bioinformatics_Assignment4_part1.Rmd   # Part 1 report (knit to HTML/PDF)
├─ Bioinformatics_Assignment4_part2.Rmd   # Part 2 report (knit to HTML/PDF)
├─ 01_part1_expression.R                  # Standalone script: Part 1 Q1–Q5
├─ 02_part1_growth.R                      # Standalone script: Part 1 Q6–Q10
├─ 03_part2_sequences.R                   # Standalone script: Part 2 Q1–Q6
├─ data/                                  # Input data (downloaded or pre-placed)
├─ figures/                               # Generated figures
├─ output/                                # Generated tables / text outputs
├─ scripts/                               # (Optional) extra code
├─ ecoli_cds.fa, cryo_cds.fa              # Large FASTA (auto-downloaded by Part 2)
└─ README.md


Note: You do not need to commit the large FASTA files (*.fa). The Part 2 Rmd will download them automatically.

Requirements

R (≥ 4.2 recommended)

R packages:

knitr, rmarkdown

seqinr

R.utils

Internet access (for downloading datasets in both parts)

The Rmds include helper code to auto-install missing packages from CRAN.

How to run
Option A — RStudio “Knit”

Open Bioinformatics_Assignment4_part1.Rmd → Knit to HTML (or PDF if you have LaTeX).

Open Bioinformatics_Assignment4_part2.Rmd → Knit.

Both reports will:

Ensure the folders data/, figures/, output/ exist.

Write outputs into figures/ and output/.

Option B — Render from Console
rmarkdown::render("Bioinformatics_Assignment4_part1.Rmd")
rmarkdown::render("Bioinformatics_Assignment4_part2.Rmd")

Inputs
Part 1

Downloads to data/:

gene_expression.tsv

growth_data.csv

Part 2

Downloads to data/, then unzips to project root (or data/ depending on your Rmd settings):

E. coli CDS FASTA (gz → fa)

Cryobacterium CDS FASTA (gz → fa)

If FASTA files already exist, they won’t be re-downloaded.

Outputs

Tables → output/

e.g., q01_first6_genes.csv, q07_mean_sd_2005_2020_by_site.csv,
p2_q5_codon_rscu.csv, p2_q6_k3_kmers_freq.csv, etc.

Figures → figures/

e.g., fig_q05_mean_hist.png, fig_q08_boxplot_2005_2020_by_site.png,
fig_p2_q3_len_boxplot.png, fig_p2_q4_nt_proportions.png,
fig_p2_q6_k3_overrepresented.png, etc.

Reports → knitted HTML/PDF in the project root.

Running the standalone scripts (optional)

If the marker wants to test sections independently (no knitting):

Part 1 – Expression (Q1–Q5)
01_part1_expression.R

Reads data/gene_expression.tsv

Writes outputs to output/ and figures/

Part 1 – Growth (Q6–Q10)
02_part1_growth.R

Reads data/growth_data.csv

Writes outputs to output/ and figures/

Part 2 – Sequences (Q1–Q6)
03_part2_sequences.R

Downloads CDS FASTA if missing

Performs counts, composition, RSCU, and k-mer analysis

Writes outputs to output/ and figures/

Run from the project root in R/RStudio:

source("01_part1_expression.R")
source("02_part1_growth.R")
source("03_part2_sequences.R")

Notes on runtime & reproducibility

Part 2 – k-mers (Q6) can take a few minutes (especially k=4,5).
If you only need a quick demo, you can temporarily edit the Rmd to run just k=3.

The code avoids re-downloading files if they already exist.

All paths are relative to the project directory for easy portability.

What to submit

Bioinformatics_Assignment4_part1.Rmd + its knitted HTML (or PDF)

Bioinformatics_Assignment4_part2.Rmd + its knitted HTML (or PDF)

Include figures/ and output/ if requested by the rubric (or provide the knitted HTMLs which embed plots).

Unless explicitly required, do not upload the large FASTA files.

Optional: .gitignore

If you push to GitHub, consider adding:

# R stuff
.Rproj.user/
.Rhistory
.RData
.Ruserdata

# Knitr/rmarkdown caches
*_cache/
cache/
*.utf8.md

# Large data
*.fa
*.fa.gz
data/*.fa
data/*.fa.gz

# Outputs
figures/
output/

Acknowledgements / References

Charif & Lobry (2007): SeqinR

Wright (1990): Effective number of codons

Deakin University Library (2025): Harvard referencing guide
