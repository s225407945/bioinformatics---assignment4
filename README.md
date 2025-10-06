
Bioinformatics Assignment 4 — Reproducible Report
==================================================

Author: Avishka Kavindi Vithanage (s225407945)  
Date: 2025-10-06

Overview
--------

This repository contains a reproducible R-based analysis for Bioinformatics Assignment 4. The assessment demonstrates the application of core bioinformatics techniques in data wrangling, statistical analysis, sequence exploration, and reproducible research practices using RMarkdown and R scripts.

The assignment is split into two parts:

- Part 1: Analysis of biological datasets including RNA-seq gene expression and tree growth measurements.
- Part 2: Comparative sequence analysis of two bacterial species — Escherichia coli and Cryobacterium sp. SO1 — using FASTA-format coding sequences.

All code is self-contained and reproducible: it automatically creates necessary folders, downloads datasets, and saves all intermediate and final outputs.

Quick Start
-----------

Run scripts individually (preferred for fast execution):

    source("01_part1_expression.R")  # Part 1 — Q1–Q5: Gene expression
    source("02_part1_growth.R")      # Part 1 — Q6–Q10: Tree growth
    source("03_part2_sequences.R")   # Part 2 — Q1–Q6: Sequence analysis

Or knit the full RMarkdown reports (for documentation + analysis):

    rmarkdown::render("Bioinformatics_Assignment4_part1.Rmd")
    rmarkdown::render("Bioinformatics_Assignment4_part2.Rmd")

Tip: Prefer HTML output. PDF works if xelatex is available and Unicode (Δ) is handled.

Folder Structure
----------------

This structure is automatically created during script or report execution:

    data/      # Input files (downloaded here)
    output/    # Output tables (CSV/TXT)
    figures/   # Plots and diagrams (PNG)
    scripts/   # (Optional) CLI helper scripts

Generic Bioinformatics Concepts
-------------------------------

Gene Expression Analysis (Part 1 — Q1–Q5)
- Process RNA-seq count data to calculate average gene expression.
- Identify highly expressed genes.
- Visualize expression distributions.

Tree Growth Analysis (Part 1 — Q6–Q10)
- Summarize growth over time.
- Visualize size distributions across sites.
- Perform statistical tests (t-tests) to compare conditions.

Sequence Diversity Analysis (Part 2 — Q1–Q6)
- Compare CDS length and GC content.
- Visualize codon usage and amino-acid composition.
- Explore codon bias and k-mer patterns.

Script Descriptions
-------------------

| File                   | Purpose                        | Inputs                 | Outputs                                       |
|------------------------|--------------------------------|------------------------|-----------------------------------------------|
| 01_part1_expression.R  | Gene expression analysis       | gene_expression.tsv    | Expression stats, histogram                   |
| 02_part1_growth.R      | Tree growth analysis           | growth_data.csv        | Growth summary, boxplots, t-tests             |
| 03_part2_sequences.R   | Sequence diversity analysis    | ecoli & cryo CDS FASTA | CDS stats, codon usage, nucleotide analysis   |

All scripts:
- Download required datasets into data/ if missing
- Auto-install required packages
- Create and populate output/ and figures/ folders

RMarkdown Report Descriptions
-----------------------------

| File                                 | Purpose                           | Outputs                                       |
|--------------------------------------|-----------------------------------|-----------------------------------------------|
| Bioinformatics_Assignment4_part1.Rmd | Narrative for Part 1 (Q1–Q10)     | All Part 1 outputs and figures                |
| Bioinformatics_Assignment4_part2.Rmd | Narrative for Part 2 (Q1–Q6)      | All Part 2 outputs and figures                |

Requirements
------------

- R ≥ 4.1
- Required R packages: rmarkdown, knitr, seqinr, R.utils
- PDF rendering (optional): LaTeX distribution with xelatex

YAML for PDF with Unicode (in .Rmd files):

    output:
      html_document: default
      pdf_document:
        latex_engine: xelatex

Reproducibility & Environment
-----------------------------

- Tested on R 4.1+ in clean environments.
- No manual paths required — everything runs from the repo root.
- Scripts are idempotent — re-running them will reuse already downloaded data and outputs.

Troubleshooting
---------------

| Issue                          | Solution                                            |
|--------------------------------|-----------------------------------------------------|
| Unicode PDF error (e.g., Δ)    | Use latex_engine: xelatex or knit to HTML          |
| No internet access             | Manually place required data files in data/        |
| Slow sequence analysis         | Only the first run is slow; uses cached data/ files|

References
----------

- Deakin University Library 2025, Deakin guide to Australian Harvard referencing.
- Charif & Lobry 2007, SeqinR 1.0-2.
- Wright 1990, “The effective number of codons used in a gene”, Gene.

