# Generate precomputed gene annotation tables for AnnotateFeatures() fallback.
#
# Run this script once (needs internet + biomaRt) to produce three .rds files.
# Upload them to the scSidekick GitHub release tagged "annotation-v1":
#
#   gh release create annotation-v1 --title "Gene annotation tables" \
#     mouse_annotation.rds human_annotation.rds zebrafish_annotation.rds
#
# Output columns (nk_* prefix, gene symbol as rowname):
#   nk_ensembl_id, nk_chr, nk_start, nk_end, nk_strand,
#   nk_biotype, nk_description

library(biomaRt)

.fetch_species <- function(dataset, label) {
  message("Fetching ", label, " from Ensembl ...")
  mart <- useMart("ensembl", dataset = dataset)

  ann <- getBM(
    attributes = c(
      "external_gene_name", "ensembl_gene_id",
      "chromosome_name",    "start_position",
      "end_position",       "strand",
      "gene_biotype",       "description"
    ),
    mart = mart
  )

  ann <- ann[nzchar(ann$external_gene_name), ]
  ann <- ann[!duplicated(ann$external_gene_name), ]
  rownames(ann) <- ann$external_gene_name
  ann$external_gene_name <- NULL

  colnames(ann) <- c(
    "nk_ensembl_id", "nk_chr", "nk_start", "nk_end",
    "nk_strand", "nk_biotype", "nk_description"
  )

  message("  → ", nrow(ann), " genes for ", label)
  ann
}

mouse     <- .fetch_species("mmusculus_gene_ensembl", "mouse")
human     <- .fetch_species("hsapiens_gene_ensembl",  "human")
zebrafish <- .fetch_species("drerio_gene_ensembl",     "zebrafish")

saveRDS(mouse,     "mouse_annotation.rds")
saveRDS(human,     "human_annotation.rds")
saveRDS(zebrafish, "zebrafish_annotation.rds")

message("Done. Upload the three .rds files to the 'annotation-v1' GitHub release.")
