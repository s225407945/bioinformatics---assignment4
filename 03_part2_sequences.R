## PART 2 — E. coli vs Cryobacterium sp. SO1
## Needs: seqinr, R.utils installed
suppressPackageStartupMessages({
  library(seqinr)
  library(R.utils)
})

outdirT <- "output"
outdirF <- "figures"
dir.create(outdirT, showWarnings = FALSE)
dir.create(outdirF, showWarnings = FALSE)

# ---- DOWNLOAD CDS FASTA (gz) ----
URL_ECOLI <- "http://ftp.ensemblgenomes.org/pub/bacteria/release-53/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/cds/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.cds.all.fa.gz"
URL_CRYO  <- "https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-62/fasta/bacteria_98_collection/cryobacterium_sp_so1_gca_004210215/cds/Cryobacterium_sp_so1_gca_004210215.ASM421021v1.cds.all.fa.gz"

if (!file.exists("ecoli_cds.fa")) {
  download.file(URL_ECOLI, destfile = "ecoli_cds.fa.gz", quiet = TRUE)
  gunzip("ecoli_cds.fa.gz", overwrite = TRUE)
}
if (!file.exists("cryo_cds.fa")) {
  download.file(URL_CRYO, destfile = "cryo_cds.fa.gz", quiet = TRUE)
  gunzip("cryo_cds.fa.gz", overwrite = TRUE)
}

# ---- READ FASTA ----
cds_ecoli <- read.fasta("ecoli_cds.fa")
cds_cryo  <- read.fasta("cryo_cds.fa")

# ---- Q1: How many CDS ----
n_ecoli <- length(cds_ecoli)
n_cryo  <- length(cds_cryo)
tab_counts <- data.frame(Organism = c("E. coli K-12 MG1655", "Cryobacterium sp. SO1"),
                         CDS_count = c(n_ecoli, n_cryo))
write.csv(tab_counts, file.path(outdirT, "p2_q1_cds_counts.csv"), row.names = FALSE)

# ---- Q2: Total coding DNA length ----
len_ecoli <- as.numeric(summary(cds_ecoli)[,1])
len_cryo  <- as.numeric(summary(cds_cryo)[,1])
tab_totbp <- data.frame(Organism = c("E. coli K-12 MG1655", "Cryobacterium sp. SO1"),
                        CDS_count = c(n_ecoli, n_cryo),
                        Total_CDS_bp = c(sum(len_ecoli), sum(len_cryo)))
write.csv(tab_totbp, file.path(outdirT, "p2_q2_total_coding_bp.csv"), row.names = FALSE)

# ---- Q3: Length distributions + boxplot + mean/median ----
len_df <- rbind(
  data.frame(Organism = "E. coli", Len = len_ecoli),
  data.frame(Organism = "Cryobacterium", Len = len_cryo)
)
stats_len <- aggregate(Len ~ Organism, data = len_df,
                       FUN = function(x) c(Mean = mean(x), Median = median(x)))
# unwrap the c() result
stats_len_out <- data.frame(Organism = stats_len$Organism,
                            Mean    = sapply(stats_len$Len, `[`, 1),
                            Median  = sapply(stats_len$Len, `[`, 2))
write.csv(stats_len_out, file.path(outdirT, "p2_q3_len_mean_median.csv"), row.names = FALSE)

png(file.path(outdirF, "fig_p2_q3_len_boxplot.png"), width = 900, height = 600)
boxplot(Len ~ Organism, data = len_df, ylab = "CDS length (bp)",
        main = "CDS length distribution")
dev.off()

# ---- Q4: Nucleotide & Amino-acid frequencies (barplots) ----
# (robust to uppercase/ambiguous bases)

# make sure FASTA are read as lists of chars (already done above)
# cds_ecoli <- read.fasta("ecoli_cds.fa", seqtype = "DNA", as.string = FALSE, forceDNAtolower = FALSE)
# cds_cryo  <- read.fasta("cryo_cds.fa",  seqtype = "DNA", as.string = FALSE, forceDNAtolower = FALSE)

# 1) Nucleotide frequencies
dna_ecoli <- tolower(unlist(cds_ecoli))
dna_cryo  <- tolower(unlist(cds_cryo))

nt_ecoli <- count(dna_ecoli, 1)
nt_cryo  <- count(dna_cryo,  1)

# keep canonical bases only (drop Ns etc.)
nt_ecoli <- nt_ecoli[c("a","c","g","t")]
nt_cryo  <- nt_cryo[c("a","c","g","t")]

nt_tab <- data.frame(
  Base = c("A","C","G","T"),
  E_coli = as.numeric(nt_ecoli),
  Cryobacterium = as.numeric(nt_cryo)
)
nt_prop <- within(nt_tab, {
  E_coli        <- E_coli / sum(E_coli)
  Cryobacterium <- Cryobacterium / sum(Cryobacterium)
})
write.csv(nt_tab,  file.path(outdirT, "p2_q4_nt_counts.csv"),      row.names = FALSE)
write.csv(nt_prop, file.path(outdirT, "p2_q4_nt_proportions.csv"), row.names = FALSE)

png(file.path(outdirF, "fig_p2_q4_nt_proportions.png"), width = 900, height = 600)
barplot(t(as.matrix(nt_prop[, -1])), beside = TRUE,
        names.arg = nt_prop$Base, ylim = c(0, max(nt_prop[, -1]) * 1.1),
        main = "Nucleotide proportions in CDS")
legend("topright", legend = c("E. coli","Cryobacterium"),
       fill = c("gray60","gray80"), bty = "n")
grid(); dev.off()

# 2) Amino-acid frequencies
aa_alphabet <- s2c("ACDEFGHIKLMNPQRSTVWY")  # 20 standard AAs
prot_ecoli <- lapply(cds_ecoli, translate)
prot_cryo  <- lapply(cds_cryo,  translate)

# concatenate and drop stop codons (*)
prot_ecoli_all <- unlist(prot_ecoli); prot_ecoli_all <- prot_ecoli_all[prot_ecoli_all != "*"]
prot_cryo_all  <- unlist(prot_cryo);  prot_cryo_all  <- prot_cryo_all[prot_cryo_all  != "*"]

aa_ecoli <- count(prot_ecoli_all, wordsize = 1, alphabet = aa_alphabet)
aa_cryo  <- count(prot_cryo_all,  wordsize = 1, alphabet = aa_alphabet)

aa_tab <- data.frame(
  AminoAcid = aa_alphabet,
  E_coli = as.numeric(aa_ecoli),
  Cryobacterium = as.numeric(aa_cryo)
)
aa_prop <- within(aa_tab, {
  E_coli        <- E_coli / sum(E_coli)
  Cryobacterium <- Cryobacterium / sum(Cryobacterium)
})
write.csv(aa_tab,  file.path(outdirT, "p2_q4_aa_counts.csv"),      row.names = FALSE)
write.csv(aa_prop, file.path(outdirT, "p2_q4_aa_proportions.csv"), row.names = FALSE)

png(file.path(outdirF, "fig_p2_q4_aa_proportions.png"), width = 1000, height = 600)
barplot(t(as.matrix(aa_prop[, -1])), beside = TRUE,
        names.arg = aa_prop$AminoAcid, las = 2,
        main = "Amino-acid proportions (translated CDS)")
legend("topright", legend = c("E. coli","Cryobacterium"),
       fill = c("gray60","gray80"), bty = "n")
grid(); dev.off()


# ---- Q5: Codon usage & bias (RSCU) ----
# Clean, concatenate, and trim to frame so length is a multiple of 3

clean_concat <- function(cds_list) {
  v <- tolower(unlist(cds_list))
  v <- v[v %in% c("a","c","g","t")]          # drop ambiguous bases
  rem <- length(v) %% 3
  if (rem != 0) v <- v[seq_len(length(v) - rem)]
  v
}
ecoli_concat <- clean_concat(cds_ecoli)
cryo_concat  <- clean_concat(cds_cryo)

rscu_ecoli <- uco(ecoli_concat, index = "rscu", as.data.frame = TRUE)
rscu_cryo  <- uco(cryo_concat,  index = "rscu", as.data.frame = TRUE)

# drop stop codons if present
rscu_ecoli <- subset(rscu_ecoli, AA != "*")
rscu_cryo  <- subset(rscu_cryo,  AA != "*")

# align codon lists
codons <- sort(unique(c(rscu_ecoli$codon, rscu_cryo$codon)))
getRSCU <- function(df, codons) {
  x <- df[match(codons, df$codon), c("codon","AA","RSCU")]
  x$RSCU[is.na(x$RSCU)] <- 0
  x
}
E <- getRSCU(rscu_ecoli, codons)
C <- getRSCU(rscu_cryo,  codons)

codon_tab <- data.frame(
  Codon = codons, AA = E$AA,
  RSCU_Ecoli = E$RSCU, RSCU_Cryo = C$RSCU,
  Diff = C$RSCU - E$RSCU
)
write.csv(codon_tab, file.path(outdirT, "p2_q5_codon_rscu.csv"), row.names = FALSE)

# plot top 20 absolute differences
ord <- order(abs(codon_tab$Diff), decreasing = TRUE)
top <- codon_tab[ord[1:20], ]
png(file.path(outdirF, "fig_p2_q5_codon_rscu_topdiff.png"), width = 1100, height = 700)
op <- par(mar = c(10,4,4,1)+0.1)
barplot(t(as.matrix(top[, c("RSCU_Ecoli","RSCU_Cryo")])),
        beside = TRUE, las = 2,
        names.arg = paste(top$Codon, paste0("(", top$AA, ")")),
        main = "Top 20 codons by |RSCU difference|")
legend("topleft", legend = c("E. coli","Cryobacterium"),
       fill = c("gray60","gray80"), bty = "n")
grid(); par(op); dev.off()


# ---- Q6: Protein k-mers (k = 3..5), over/under in Cryo vs E. coli ----
# Add a light check + plots with grid()
stopifnot(length(prot_ecoli_all) > 0, length(prot_cryo_all) > 0, length(aa_alphabet) == 20)

compare_kmers <- function(k) {
  k_ec <- count(prot_ecoli_all, wordsize = k, alphabet = aa_alphabet, freq = TRUE)
  k_cr <- count(prot_cryo_all,  wordsize = k, alphabet = aa_alphabet, freq = TRUE)
  all_names <- sort(unique(c(names(k_ec), names(k_cr))))
  k1 <- k_ec[all_names]; k1[is.na(k1)] <- 0
  k2 <- k_cr[all_names]; k2[is.na(k2)] <- 0
  
  df <- data.frame(Kmer = all_names,
                   Ecoli = as.numeric(k1),
                   Cryo  = as.numeric(k2))
  df$Diff <- df$Cryo - df$Ecoli
  df <- df[order(df$Diff, decreasing = TRUE), ]
  
  write.csv(df, file.path(outdirT, sprintf("p2_q6_k%d_kmers_freq.csv", k)), row.names = FALSE)
  
  top10 <- head(df, 10)
  bot10 <- tail(df, 10)
  
  png(file.path(outdirF, sprintf("fig_p2_q6_k%d_overrepresented.png", k)),
      width = 1100, height = 700)
  barplot(top10$Diff, names.arg = top10$Kmer, las = 2,
          main = sprintf("k=%d: Top 10 over-represented in Cryobacterium (Cryo − E. coli)", k),
          ylab = "Frequency difference")
  grid(); dev.off()
  
  png(file.path(outdirF, sprintf("fig_p2_q6_k%d_underrepresented.png", k)),
      width = 1100, height = 700)
  barplot(rev(bot10$Diff), names.arg = rev(bot10$Kmer), las = 2,
          main = sprintf("k=%d: Top 10 under-represented in Cryobacterium (Cryo − E. coli)", k),
          ylab = "Frequency difference")
  grid(); dev.off()
}

for (k in 3:5) compare_kmers(k)
