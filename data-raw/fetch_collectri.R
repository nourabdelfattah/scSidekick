# data-raw/fetch_collectri.R
# Run this script ONCE (when OmniPath is reachable) to bundle CollecTRI
# as inst/extdata/collectri_{human,mouse}.rds so RunURA() never needs the
# OmniPath server at runtime.
#
# Usage:
#   source("data-raw/fetch_collectri.R")

.fetch_one <- function(organism, taxon) {
  message("Fetching CollecTRI for ", organism, " (taxon ", taxon, ")...")

  # Attempt 1: decoupleR wrapper
  net <- tryCatch(
    suppressMessages(decoupleR::get_collectri(organism = organism,
                                              split_complexes = FALSE)),
    error = function(e) NULL
  )

  # Attempt 2: direct OmniPath REST API (bypasses OmnipathR entirely)
  if (is.null(net) || nrow(net) == 0L) {
    message("  decoupleR/OmnipathR failed — trying direct REST download...")
    url <- paste0(
      "https://omnipathdb.org/interactions",
      "?datasets=collectri",
      "&organisms=", taxon,
      "&genesymbols=1",
      "&fields=source_genesymbol,target_genesymbol,is_stimulation,is_inhibition"
    )
    tmp <- tempfile(fileext = ".tsv")
    ok  <- tryCatch({
      utils::download.file(url, tmp, quiet = TRUE)
      TRUE
    }, error = function(e) FALSE)

    if (ok) {
      raw <- utils::read.table(tmp, header = TRUE, sep = "\t",
                               quote = "", stringsAsFactors = FALSE)
      unlink(tmp)
      src <- intersect(c("source_genesymbol", "source"), colnames(raw))[1L]
      tgt <- intersect(c("target_genesymbol", "target"), colnames(raw))[1L]
      if (!is.na(src) && !is.na(tgt)) {
        mor <- if ("mor" %in% colnames(raw)) {
          as.numeric(raw$mor)
        } else if (all(c("is_stimulation", "is_inhibition") %in% colnames(raw))) {
          ifelse(as.integer(raw$is_stimulation) == 1L, 1, -1)
        } else {
          rep(1, nrow(raw))
        }
        net <- data.frame(source = as.character(raw[[src]]),
                          target = as.character(raw[[tgt]]),
                          mor    = mor,
                          stringsAsFactors = FALSE)
        net <- net[nzchar(net$source) & nzchar(net$target) & is.finite(net$mor), ]
      }
    }
  }

  if (is.null(net) || nrow(net) == 0L)
    stop("Could not fetch CollecTRI for ", organism,
         ". Check network access or try again later.")

  n_tfs <- length(unique(net$source))
  message("  OK — ", nrow(net), " interactions, ", n_tfs, " TFs.")
  net[, c("source", "target", "mor")]
}

# Fetch and save
for (spec in list(list(org = "human", taxon = "9606"),
                  list(org = "mouse", taxon = "10090"))) {
  net  <- .fetch_one(spec$org, spec$taxon)
  path <- file.path("inst", "extdata",
                    paste0("collectri_", spec$org, ".rds"))
  saveRDS(net, path)
  message("Saved → ", path)
}

message("\nDone. Re-install scSidekick (devtools::load_all()) to pick up the bundled data.")
