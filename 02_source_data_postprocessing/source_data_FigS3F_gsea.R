# PAPER PANEL: Fig S3F (source data)
# PLOTS:  the 12 GO biological-process pathways shown in Fig S3F, chosen from the significant GSEA results
#         (all with adjusted p < 0.001): 7 higher in aged P1KO (rescue) and 5 higher in aged CTRL (aging)
# CELLS:  Aged CTRL and Aged P1KO visceral SMCs
# GROUPS COMPARED: Aged P1KO vs Aged CTRL
# TEST:   GSEA (fgsea) on the genome-wide ranking; adjusted p
library(openxlsx)
xlsx <- "output/source_data/FigS3F_data.xlsx"
all_sig <- read.xlsx(xlsx, sheet = "all_significant")
de_sum  <- read.xlsx(xlsx, sheet = "de_summary")
shown <- c("GOBP_SARCOPLASMIC_RETICULUM_CALCIUM_ION_TRANSPORT", "GOBP_CARDIAC_MUSCLE_CONTRACTION",
           "GOBP_STRIATED_MUSCLE_CONTRACTION", "GOBP_MUSCLE_CONTRACTION", "GOBP_MUSCLE_SYSTEM_PROCESS",
           "GOBP_HEART_PROCESS", "GOBP_MUSCLE_STRUCTURE_DEVELOPMENT",
           "GOBP_DEFENSE_RESPONSE_TO_OTHER_ORGANISM", "GOBP_IMMUNE_EFFECTOR_PROCESS",
           "GOBP_DEFENSE_RESPONSE_TO_SYMBIONT", "GOBP_RIBOSOME_BIOGENESIS", "GOBP_HUMORAL_IMMUNE_RESPONSE")
stopifnot(all(shown %in% all_sig$pathway))
d <- all_sig[match(shown, all_sig$pathway), c("pathway", "NES", "padj", "size")]
d$pathway_clean <- tools::toTitleCase(tolower(gsub("_", " ", sub("^GOBP_", "", d$pathway))))
d$direction <- ifelse(d$NES > 0, "Higher in KO (rescue)", "Higher in Ctrl (aging)")
d <- d[, c("pathway", "pathway_clean", "NES", "padj", "size", "direction")]
wb <- createWorkbook()
addWorksheet(wb, "go_enrichment"); writeData(wb, "go_enrichment", d)
addWorksheet(wb, "all_significant"); writeData(wb, "all_significant", all_sig)
addWorksheet(wb, "de_summary"); writeData(wb, "de_summary", de_sum)
saveWorkbook(wb, xlsx, overwrite = TRUE)
cat("Fig S3F source data:", nrow(d), "pathways\n")
