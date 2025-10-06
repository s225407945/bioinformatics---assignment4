## PART 1 — gene_expression
## Reads counts, computes stats, saves tables & plots (base R only)

# paths
infile  <- "data/gene_expression.tsv"
outdirT <- "output"
outdirF <- "figures"

# read (gene IDs as rownames)
gene <- read.table(infile, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# Q1: first six genes as a small table
q1_head <- head(gene)
write.csv(q1_head, file.path(outdirT, "q01_first6_genes.csv"))

# Add mean column
gene$mean_expr <- rowMeans(gene)

# Q2: first six with mean
write.csv(head(gene), file.path(outdirT, "q02_first6_with_mean.csv"))

# Q3: top 10 by mean
o <- order(gene$mean_expr, decreasing = TRUE)
top10 <- head(gene[o, , drop = FALSE], 10)
write.csv(top10, file.path(outdirT, "q03_top10_highest_mean.csv"))

# Q4: number of genes with mean < 10
n_lt10 <- sum(gene$mean_expr < 10)
writeLines(sprintf("Genes with mean < 10: %d", n_lt10),
           con = file.path(outdirT, "q04_count_mean_lt10.txt"))

# Q5: histogram plot
png(file.path(outdirF, "fig_q05_mean_hist.png"), width = 900, height = 600)
hist(gene$mean_expr, breaks = 40,
     xlab = "Mean expression", main = "Q5: Histogram of mean expression")
grid()
dev.off()
