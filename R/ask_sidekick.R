# =============================================================================
# scSidekick - AskSidekick()  (ask_sidekick.R)
#
# A lightweight AI assistant for navigating scSidekick: "which function do I
# need", "which metadata column looks like X", "how do I recreate this
# figure". Backend: Google Gemini's free tier via the 'ellmer' package
# (Suggests only - never required just to load scSidekick).
#
# Each call is independent (no multi-turn memory): pasting the package's
# entire documentation into every call would blow through the free tier's
# per-minute token limit in ~1-2 calls, so instead each question retrieves
# only its relevant documentation chunks (keyword overlap + exact
# function-name match).
#
# The documentation context is built from the INSTALLED package itself
# (tools::Rd_db() for man pages, system.file("doc") for vignette sources -
# both work against any real install, not just a source checkout) and cached
# under tools::R_user_dir("scSidekick", "cache"), keyed by package version.
# The first AskSidekick() call after installing/updating scSidekick rebuilds
# it automatically; every call after that reuses the cache. No bundled
# asset, no manual release-day rebuild step - RebuildSidekickContext() is
# there if you want to force it (e.g. while actively developing scSidekick
# itself, where docs can change without a version bump).
#
# Free-tier daily quotas (RPD) vary a lot by model - "Flash Lite" models
# tend to get a much larger daily allowance than flagship Flash models on
# the same account (check https://aistudio.google.com/rate-limit for yours).
# Each model is a separate quota pool, so instead of waiting out a 429,
# AskSidekick() tries the next model in `models` and says so when it does
# (never a silent downgrade).
# =============================================================================

.sidekick_context_version <- function() as.character(utils::packageVersion("scSidekick"))

.sidekick_cache_path <- function() file.path(.nk_cache_dir(), "sidekick_context.rds")

.sidekick_default_models <- c(
  "gemini-3.5-flash-lite",
  "gemini-3.1-flash-lite",
  "gemini-3.7-flash",
  "gemini-3.5-flash",
  "gemini-3-flash-preview",
  "gemini-2.5-flash"
)

.sidekick_system_preamble <- paste0(
  "You are Sidekick, a help assistant for the scSidekick R package (a ",
  "single-cell / spatial transcriptomics toolkit built on Seurat). Your ",
  "users are often new to R.\n\n",
  "Rules:\n",
  "- Only reference functions, arguments, and defaults that actually appear ",
  "in the documentation excerpts provided below. Never invent a function ",
  "name or argument. If the excerpts don't cover something, say so ",
  "explicitly instead of guessing.\n",
  "- Prefer giving a short, copy-pasteable R code snippet over a long ",
  "explanation.\n",
  "- If a specific variable name for the user's object is given below, use ",
  "that exact name in every code snippet - never substitute a generic ",
  "placeholder like `seurat_obj` or `obj`.\n",
  "- Call every scSidekick function as `scSidekick::FunctionName(...)` in code ",
  "snippets, rather than assuming `library(scSidekick)` has already been run. ",
  "This makes each snippet self-contained and unambiguous about which package ",
  "a function comes from.\n",
  "- Keep answers brief - a sentence or two of explanation plus code, not an ",
  "essay.\n",
  "- Several scSidekick plotting functions only work on a column or data ",
  "structure first produced by a prerequisite `Run*()` function - e.g. ",
  "`PlotPseudotime()`/`PlotFeatureTrend()`/`PlotTrajectory()` need a ",
  "pseudotime column from `RunSlingshot()`; pathway `Plot*()` functions need ",
  "`RunGSEA()`/`RunGSEA_pseudobulk()`/`RunSCssGSEA()` to have already run. ",
  "If the user's object doesn't yet have the column/data a function needs ",
  "(check the documentation excerpts for what each function actually reads), ",
  "never substitute an existing but wrong-typed column (e.g. a categorical ",
  "column standing in for a numeric pseudotime value) just to fill in the ",
  "call. Instead, ADD THE PREREQUISITE CALL AS AN ACTUAL LINE INSIDE THE SAME ",
  "code block, immediately before the line that needs its output - not as a ",
  "side note in the surrounding prose, where it's easy to miss. The snippet ",
  "should read top to bottom as a runnable sequence: prerequisite call(s) ",
  "first, assigning back onto the user's own object variable, then the ",
  "plotting call that consumes what it produced.\n"
)

.sidekick_stopwords <- c(
  "the","a","an","is","are","do","does","how","what","which","use","using",
  "used","function","for","to","of","in","on","with","and","or","i","my",
  "this","that","it","can","would","should","need","want","get","find","from"
)

# ── Documentation context: build from the installed package, cache by version ──

# Man pages via tools::Rd_db() (works against any installed copy, not just a
# source checkout) and vignette SOURCES via system.file("doc") (requires
# VignetteBuilder: knitr in DESCRIPTION, which builds inst/doc/*.Rmd on
# install). Rd2txt() renders bold/italic as terminal-style overstrike
# ("_\bS_\bt..."); stripped below since it roughly doubles token count for no
# benefit here.
.sidekick_build_context_from_install <- function() {
  chunks <- list()
  add_chunk <- function(id, label, text, always = FALSE) {
    text <- gsub(".\b", "", text, perl = TRUE)
    # Force an explicit UTF-8 encoding mark regardless of the machine's
    # locale (vignettes/.Rd text can contain non-ASCII characters like en
    # dashes, arrows, "x", "<="). Without this, text read under a non-UTF-8
    # locale (e.g. the "C" locale) carries no/incorrect encoding metadata,
    # which surfaces later as "input string is invalid in this locale"
    # warnings out of readRDS() or string ops on the cached bundle.
    text <- enc2utf8(text)
    if (nzchar(trimws(text)))
      chunks[[id]] <<- list(label = label, text = text, always = always)
  }

  # Extract just a function's one-line \title{} from its parsed Rd object -
  # avoids Rd2txt (which renders the whole page) for the tiny always-included
  # index chunk below. "title" is a safe, escaping-free substring match for
  # the Rd_tag (no other tag - \name, \alias, \description, \usage,
  # \arguments, \value - contains it).
  .rd_title <- function(rd) {
    parts <- Filter(function(x) grepl("title", attr(x, "Rd_tag") %||% "",
                                      fixed = TRUE), rd)
    if (!length(parts)) return("")
    trimws(paste(vapply(parts, function(p) paste(as.character(p), collapse = ""),
                        character(1)), collapse = " "))
  }

  db <- tryCatch(tools::Rd_db("scSidekick"), error = function(e) list())

  # ── Always-included function index ──────────────────────────────────────
  # A one-line-per-function index (name + \title), ALWAYS placed in the
  # system prompt regardless of what keyword retrieval scores highest.
  # Without this, a question whose best-matching man page retrieval misses
  # (imperfect keyword overlap, or a path-matching question whose own
  # instructional text happens to name OTHER functions as illustrative
  # examples - see path_guidance below) leaves the model with literally no
  # signal about which functions even exist, and it either invents one or
  # falls back to generic Seurat/ggplot code instead of a real scSidekick
  # function. Cheap (under ~5KB) and single-sourced from the same roxygen
  # titles already used for the man pages, so it can't drift out of sync.
  index_lines <- vapply(names(db), function(nm) {
    fn_name <- sub("\\.Rd$", "", nm)
    title <- tryCatch(.rd_title(db[[nm]]), error = function(e) "")
    if (!nzchar(title)) return("")
    paste0("- `", fn_name, "()`: ", title)
  }, character(1))
  index_lines <- index_lines[nzchar(index_lines)]
  if (length(index_lines) > 0) {
    add_chunk("function_index", "scSidekick function index",
              paste(sort(index_lines), collapse = "\n"), always = TRUE)
  }

  for (nm in names(db)) {
    fn_name <- sub("\\.Rd$", "", nm)
    txt <- tryCatch({
      # Close explicitly right here rather than via on.exit(): tryCatch({...})
      # evaluates in THIS function's own frame (not a fresh call frame per
      # iteration), so on.exit() inside a loop body would replace itself each
      # time (no add = TRUE) and never actually close a connection until the
      # whole function returns - leaving every later iteration's
      # textConnection("out", ...) call colliding with a still-open one from
      # the iteration before it.
      con <- textConnection("out", "w", local = TRUE)
      tools::Rd2txt(db[[nm]], out = con, package = "scSidekick")
      close(con)
      paste(out, collapse = "\n")
    }, error = function(e) "")
    add_chunk(paste0("man:", fn_name), paste0("man page: ", fn_name), txt)
  }

  doc_dir <- system.file("doc", package = "scSidekick")
  if (nzchar(doc_dir)) {
    for (f in list.files(doc_dir, pattern = "\\.Rmd$", full.names = TRUE)) {
      # Open an explicit UTF-8 connection (.Rmd sources are always saved as
      # UTF-8) rather than a bare readLines(), which would interpret the raw
      # bytes using whatever the reading machine's locale defaults to.
      con <- file(f, encoding = "UTF-8")
      txt <- tryCatch(paste(readLines(con, warn = FALSE), collapse = "\n"),
                      error = function(e) "")
      close(con)
      add_chunk(paste0("vignette:", basename(f)), paste0("vignette: ", basename(f)), txt)
    }
  }

  chunks
}

.sidekick_get_context <- function(rebuild = FALSE) {
  cache_path <- .sidekick_cache_path()
  current_version <- .sidekick_context_version()

  if (!rebuild && file.exists(cache_path)) {
    cached <- tryCatch(readRDS(cache_path), error = function(e) NULL)
    if (!is.null(cached) && identical(cached$version, current_version) &&
        length(cached$chunks) > 0) {
      return(cached$chunks)
    }
  }

  message("scSidekick: setting up AskSidekick() for the first time with this ",
          "package version (v", current_version, ") - hang on, this can take ",
          "a minute or two. It only happens once; every call after this is ",
          "instant.")
  chunks <- .sidekick_build_context_from_install()
  dir.create(.nk_cache_dir(), recursive = TRUE, showWarnings = FALSE)
  tryCatch(
    saveRDS(list(version = current_version, built_at = Sys.time(), chunks = chunks),
           cache_path),
    error = function(e)
      warning("Could not cache AskSidekick() context: ", conditionMessage(e))
  )
  message("scSidekick: ready (", length(chunks), " documentation chunks cached).")
  chunks
}

#' Force a rebuild of AskSidekick()'s documentation cache
#'
#' `AskSidekick()` automatically rebuilds its cached documentation context
#' the first time it's called after a scSidekick version change, so this is
#' rarely needed. Use it if you're actively developing scSidekick itself
#' (docs changed without a version bump) or just want to be sure the cache
#' reflects the currently installed package.
#'
#' @return Invisibly, the number of documentation chunks cached.
#' @export
RebuildSidekickContext <- function() {
  chunks <- .sidekick_get_context(rebuild = TRUE)
  message("scSidekick: AskSidekick() context rebuilt (", length(chunks),
          " documentation chunks).")
  invisible(length(chunks))
}

# ── Retrieval: exact function-name mentions win big; otherwise score chunks ──
# by how many distinct question keywords appear in them. Cheap, local, no
# extra API calls, and good enough when questions usually name a task or a
# function.
.sidekick_retrieve <- function(question, chunks, extra_keywords = character(0),
                               token_budget = 20000L) {
  words <- unique(tolower(unlist(strsplit(question %||% "", "[^A-Za-z0-9_.]+"))))
  words <- words[nchar(words) >= 4 & !words %in% .sidekick_stopwords]
  words <- c(words, tolower(extra_keywords))

  fn_mentions <- unique(unlist(regmatches(
    question %||% "", gregexpr("[A-Z][A-Za-z0-9_.]{2,}", question %||% "")
  )))

  scored <- vapply(names(chunks), function(id) {
    ch <- chunks[[id]]
    if (isTRUE(ch$always)) return(-1)  # always-chunks are handled separately
    score <- 0
    if (startsWith(id, "man:")) {
      fn_name <- sub("^man:", "", id)
      if (fn_name %in% fn_mentions) score <- score + 100
    }
    txt_lower <- tolower(ch$text)
    score <- score + sum(vapply(words, function(w) grepl(w, txt_lower, fixed = TRUE),
                                logical(1)))
    score
  }, numeric(1))

  ordered <- names(sort(scored[scored > 0], decreasing = TRUE))

  budget_chars <- token_budget * 4L
  used <- 0L
  keep <- character(0)
  for (id in ordered) {
    n <- nchar(chunks[[id]]$text)
    if (used + n > budget_chars && length(keep) > 0) break
    keep <- c(keep, id)
    used <- used + n
  }
  keep
}

.sidekick_build_system_prompt <- function(question, chunks, extra_keywords = character(0)) {
  # Always-included chunks (currently just the function index) ride along
  # regardless of retrieval score - see .sidekick_build_context_from_install().
  always_ids   <- names(chunks)[vapply(chunks, function(x) isTRUE(x$always), logical(1))]
  retrieved_ids <- .sidekick_retrieve(question, chunks, extra_keywords)
  use_ids      <- unique(c(always_ids, retrieved_ids))

  docs_text <- if (length(use_ids) == 0) {
    "(No specific documentation matched this question.)"
  } else {
    paste(vapply(use_ids, function(id) {
      paste0("\n\n----- ", chunks[[id]]$label, " -----\n\n", chunks[[id]]$text)
    }, character(1)), collapse = "")
  }

  paste0(.sidekick_system_preamble,
         "\n\n===== scSidekick documentation (relevant excerpts) =====",
         docs_text)
}

# ── Model fallback ───────────────────────────────────────────────────────────
# Remembers, within this R session, which model index last worked - so once a
# model is known to be over quota for the day, later calls skip straight past
# it instead of spending a request rediscovering that every time.
.sidekick_fallback_env <- new.env(parent = emptyenv())
.sidekick_fallback_env$start_idx <- 1L

# Cycle through `models` on any failure (quota or transient server error)
# rather than waiting: separate models have separate quota pools, so
# switching is both faster and more reliable than a blind Sys.sleep, and
# every model change is printed so it's never a silent, confusing downgrade.
.sidekick_chat_with_fallback <- function(system_prompt, parts, key, models) {
  n <- length(models)
  start <- min(.sidekick_fallback_env$start_idx, n)

  last_error <- NULL
  for (i in seq(start, n)) {
    model <- models[i]
    # echo = "none": ellmer's own default (echo = NULL) decides whether to
    # print based on whether the IMMEDIATE calling frame looks like a
    # top-level user console call - which it never is here, since chat$chat()
    # is invoked from inside this package's own wrapper function, not
    # directly from the user's console. That heuristic silently suppresses
    # all output no matter how AskSidekick() itself is called, which is
    # exactly backwards for a function whose entire purpose is showing the
    # user an answer. Disable it entirely and print explicitly, below,
    # instead - guaranteed and consistent regardless of call-stack depth.
    chat <- ellmer::chat_google_gemini(system_prompt = system_prompt,
                                       model = model, api_key = key,
                                       echo = "none")
    result <- tryCatch(list(ok = TRUE, value = do.call(chat$chat, parts)),
                       error = function(e) list(ok = FALSE, error = e))
    if (isTRUE(result$ok)) {
      if (i > start) message("scSidekick bot: answered using ", model,
                             " (earlier model(s) were over today's quota).")
      .sidekick_fallback_env$start_idx <- i
      return(result$value)
    }
    message("scSidekick bot: ", model, " unavailable right now (",
            sub("\\n.*", "", conditionMessage(result$error)), "), ",
            if (i < n) "trying the next model..." else "no models left.")
    last_error <- result$error
  }

  stop("scSidekick bot: every configured model is over quota or unavailable ",
       "today. Wait until tomorrow (daily quotas reset), or pass a different ",
       "`models` list.\n(Last error: ", conditionMessage(last_error), ")")
}

# Classifies a vector of gene symbols as "human" (ALL-CAPS), "mouse"
# (Title-case), or "mixed"/"" (no clear convention). Used only in the
# post-hoc .sidekick_check_gene_matches() check below - genes read off an
# image are always transcribed literally, as seen, never pre-converted to
# match the object's own convention, so what's reported never silently
# diverges from what's actually visible in the figure.
.sidekick_gene_case_style <- function(genes) {
  genes <- genes[nchar(genes) > 1]
  if (length(genes) == 0) return("")
  upper_frac <- mean(genes == toupper(genes))
  if (upper_frac > 0.7) "human" else if (upper_frac < 0.3) "mouse" else "mixed"
}

# ── Seurat object introspection ─────────────────────────────────────────────
# Summarizes meta.data (name, type, a few example values) and any colors/
# settings already stored by PrepObject(), so the bot can answer "which
# column is Age?" / "what colors do I already have saved?" against the
# TRAINEE'S ACTUAL OBJECT rather than guessing from column names alone.
.sidekick_describe_object <- function(seurat_object, var_name = "seurat_object") {
  meta <- seurat_object@meta.data
  col_lines <- vapply(colnames(meta), function(cn) {
    v <- meta[[cn]]
    type <- if (is.numeric(v)) "numeric" else if (is.factor(v)) "factor" else "character"
    examples <- tryCatch({
      u <- if (is.factor(v)) levels(v) else unique(stats::na.omit(v))
      paste(utils::head(u, 5), collapse = ", ")
    }, error = function(e) "")
    sprintf("  - %s (%s): %s", cn, type, examples)
  }, character(1))

  stored_colors <- tryCatch(seurat_object@misc$nk_settings$colors, error = function(e) NULL)
  color_lines <- if (length(stored_colors) > 0) {
    paste0("Variables with colors already saved via PrepObject(): ",
           paste(names(stored_colors), collapse = ", "))
  } else {
    "No colors have been saved via PrepObject() on this object yet."
  }

  paste0(
    "The user's Seurat object is assigned to the R variable named `", var_name,
    "` in their own session. ALWAYS use exactly that variable name (`", var_name,
    "`) in any R code you write - never substitute a generic placeholder like ",
    "`seurat_obj`, `obj`, or `seurat_object`.\n\n",
    "That object's meta.data columns are:\n",
    paste(col_lines, collapse = "\n"),
    "\n\n", color_lines
  )
}

# ── Gene-name sanity check ──────────────────────────────────────────────────
# The bot sometimes proposes gene names (e.g. read off a heatmap image) that
# don't actually exist in the user's object - most often because the image
# was made from data of a different species (mouse "Stat3" vs human "STAT3")
# but also just typos or genes not in this assay. Rather than have the bot
# (or a case-convention heuristic) GUESS at species, check the one thing that
# actually matters: do the proposed genes exist in rownames(seurat_object)?
# That's a direct, deterministic answer, not a guess - and it generalizes
# beyond species mismatches to any reason the suggested code wouldn't work.
.sidekick_check_gene_matches <- function(reply, seurat_object, var_name) {
  if (is.null(seurat_object)) return(NULL)

  # Pull quoted tokens ONLY from genuine c("A", "B", ...) vector literals (2+
  # comma-separated quoted items) - NOT any quoted string anywhere in the
  # reply. Single quoted values used as ordinary argument values (e.g.
  # group.by = "Subclass") are real column names or paths, not gene lists,
  # and would otherwise falsely trigger this check on every reply that names
  # a metadata column in code.
  vec_literals <- regmatches(reply, gregexpr(
    'c\\(\\s*"[^"]+"\\s*(,\\s*"[^"]+"\\s*)+\\)', reply))[[1]]
  candidates <- unique(unlist(regmatches(vec_literals, gregexpr('"[^"]+"', vec_literals))))
  candidates <- gsub('"', "", candidates)
  candidates <- candidates[grepl("^[A-Za-z][A-Za-z0-9._-]*$", candidates)]
  if (length(candidates) < 3) return(NULL)  # not enough to say anything useful

  obj_genes <- rownames(seurat_object)
  matched <- candidates %in% obj_genes
  match_rate <- mean(matched)
  if (match_rate >= 0.5) return(NULL)  # good enough match, no warning needed

  sample_obj_genes <- utils::head(obj_genes, 200)
  obj_style  <- .sidekick_gene_case_style(sample_obj_genes)
  cand_style <- .sidekick_gene_case_style(candidates)
  example    <- sample_obj_genes[nchar(sample_obj_genes) > 1][1]
  vec_txt    <- paste0('c(', paste0('"', candidates, '"', collapse = ', '), ')')

  # A case difference (mouse "Stat3" vs. human "STAT3") is the visible
  # symptom, not itself a safe fix: some orthologs share nothing beyond case,
  # others have unrelated names entirely, so this package NEVER guesses a
  # "corrected" list on its own. When the mismatch looks species-related, it
  # instead hands back a real, runnable conversion using babelgene - a
  # lightweight, offline ortholog table (no network call, no bundled
  # dependency added here; it's only ever suggested as text) - so the user
  # gets a verifiable answer instead of a silent guess.
  species_note <- if (obj_style == "human" && cand_style == "mouse") {
    paste0(
      " Your object's gene names look ALL-CAPS (human convention: e.g. \"",
      example, "\"), while the suggested genes look Title-case (mouse ",
      "convention) - this looks like a species mismatch, not a typo. ",
      "scSidekick doesn't guess at ortholog conversions itself (case alone ",
      "isn't a reliable rule - some orthologs differ by more than case), but ",
      "here's a script to convert these mouse symbols to their human ",
      "orthologs via the `babelgene` package (offline, no bundled ",
      "dependency - install once if needed):\n\n```r\n",
      "# install.packages(\"babelgene\")\n",
      "mouse_genes <- ", vec_txt, "\n",
      "converted   <- babelgene::orthologs(genes = mouse_genes, species = \"mouse\", human = FALSE)\n",
      "human_genes <- converted$human_symbol\n",
      "human_genes  # use this vector in place of the mouse symbols above\n",
      "# genes with no listed ortholog: setdiff(mouse_genes, converted$symbol)\n",
      "```\n"
    )
  } else if (obj_style == "mouse" && cand_style == "human") {
    paste0(
      " Your object's gene names look Title-case (mouse convention: e.g. \"",
      example, "\"), while the suggested genes look ALL-CAPS (human ",
      "convention) - this looks like a species mismatch, not a typo. ",
      "scSidekick doesn't guess at ortholog conversions itself (case alone ",
      "isn't a reliable rule - some orthologs differ by more than case), but ",
      "here's a script to convert these human symbols to their mouse ",
      "orthologs via the `babelgene` package (offline, no bundled ",
      "dependency - install once if needed):\n\n```r\n",
      "# install.packages(\"babelgene\")\n",
      "human_genes <- ", vec_txt, "\n",
      "converted   <- babelgene::orthologs(genes = human_genes, species = \"mouse\", human = TRUE)\n",
      "mouse_genes <- converted$symbol\n",
      "mouse_genes  # use this vector in place of the human symbols above\n",
      "# genes with no listed ortholog: setdiff(human_genes, converted$human_symbol)\n",
      "```\n"
    )
  } else paste0(
    " scSidekick doesn't currently have a built-in ortholog-conversion ",
    "function, so double-check these against `rownames(", var_name,
    ")` (or a resource like `babelgene`/biomaRt) before running the code above."
  )

  paste0(
    "\n\n---\n**Heads up**: only ", sum(matched), " of ", length(candidates),
    " suggested gene names (", round(match_rate * 100), "%) were found in `",
    var_name, "`'s actual gene names, exactly as transcribed from the image.",
    species_note
  )
}

# Generalizes the same idea to metadata COLUMN names: the bot often names a
# column in a `group.by = "X"` / `stat.by = "X"` / `fill_variable = "X"`
# style argument - check whether "X" literally exists in the user's object.
# Deterministic, same as the gene check; `same_dataset` only changes the TONE
# of the message (a miss is a real discrepancy worth flagging hard when the
# user confirmed it's the same object, vs. an expected rough edge when it's
# genuinely from elsewhere).
.sidekick_check_column_matches <- function(reply, seurat_object, var_name, same_dataset) {
  if (is.null(seurat_object)) return(NULL)

  col_args <- regmatches(reply, gregexpr(
    '[A-Za-z_.]*(\\.by|[Vv]ariables?|[Cc]olumns?)\\s*=\\s*"[^"]+"', reply))[[1]]
  if (length(col_args) == 0) return(NULL)
  candidates <- unique(sub('.*=\\s*"([^"]+)"', "\\1", col_args))
  candidates <- candidates[nzchar(candidates)]
  if (length(candidates) == 0) return(NULL)

  real_cols <- colnames(seurat_object@meta.data)
  missing <- candidates[!candidates %in% real_cols]
  if (length(missing) == 0) return(NULL)  # everything checks out

  closest_for <- function(x) {
    d <- utils::adist(tolower(x), tolower(real_cols))[1, ]
    paste(real_cols[order(d)][seq_len(min(3, length(real_cols)))], collapse = ", ")
  }
  detail <- paste(sprintf('  - "%s" not found. Closest real column(s): %s',
                          missing, vapply(missing, closest_for, character(1))),
                  collapse = "\n")

  if (isTRUE(same_dataset)) {
    paste0(
      "\n\n---\n**Heads up**: you said this figure was made from `", var_name,
      "` itself, so an exact column match should exist for every name used, ",
      "but it doesn't for:\n", detail,
      "\nEither the wrong column was picked, or it's been renamed/removed ",
      "since the figure was made - double-check before running the code above."
    )
  } else {
    paste0(
      "\n\n---\n**Note**: since this figure's source dataset wasn't confirmed ",
      "to be `", var_name, "` itself, these column names are the bot's best ",
      "estimate, not guaranteed matches:\n", detail,
      "\nVerify against the actual figure before running the code above."
    )
  }
}

# ── Main entry point ─────────────────────────────────────────────────────────
#' Ask the scSidekick help assistant a question
#'
#' A lightweight AI assistant grounded in scSidekick's own documentation
#' (rebuilt automatically from your installed copy, so it never quotes a
#' stale function signature). Answers "which function do I need", identifies
#' likely metadata columns from your actual object, and can read a heatmap
#' image to extract gene names or work backward from a saved figure to the
#' function/arguments that likely made it.
#'
#' Requires a free Google Gemini API key: get one at
#' \url{https://aistudio.google.com/apikey}, then add a line to
#' \code{~/.Renviron}: \code{GEMINI_API_KEY=your_key_here} and restart R.
#'
#' Each call is independent (no memory of earlier questions) - the free
#' tier's token-per-minute limit gets used up in 1-2 calls if full
#' conversation history were kept, so ask one self-contained question at a
#' time.
#'
#' @param question Character. What you want to ask. Optional when
#'   `image_path` is supplied (a sensible default prompt is used).
#' @param seurat_object A Seurat object, optional. When supplied, its
#'   meta.data columns and any `PrepObject()`-saved colors are given to the
#'   assistant as extra context (e.g. for "which column is Age?" or "what
#'   colors do I already have?"), and any column/gene names it proposes are
#'   checked against the object's real names afterward.
#' @param image_path Character path to an image (PNG/JPEG or PDF), optional.
#'   Use this for a screenshot or saved PDF of a figure: the assistant reads
#'   off row/column labels, identifies the likely scSidekick function (using
#'   the file's path - scSidekick functions often name their own output
#'   folders/files distinctively, e.g. `RunSCssGSEA()`'s `"SC ssGSEA <label>"`
#'   subfolder), and proposes the call to recreate it.
#' @param same_dataset Logical, `NA` by default. Only matters when
#'   `image_path` + `seurat_object` are both given: `TRUE` means the figure
#'   was made from THIS exact object, so the assistant should demand an exact
#'   metadata-column match rather than settling for an approximate one (many
#'   objects have several similar-but-not-identical columns, e.g. `Cluster`
#'   vs. `seurat_clusters` vs. `Assignment`). `FALSE`/`NA` means the figure
#'   may be from elsewhere (a different object, an earlier package version, a
#'   collaborator's file), so the assistant gives its best estimate and says
#'   so explicitly rather than implying a confirmed match.
#' @param models Character vector of Gemini model names to try, in order.
#'   Defaults to a mix of "Flash Lite" and flagship Flash models, since free
#'   tier daily quotas (RPD) differ a lot by model and each is a separate
#'   quota pool - see `https://aistudio.google.com/rate-limit` for your
#'   account's actual per-model limits, which vary by account. On a 429/503,
#'   `AskSidekick()` automatically tries the next model in this list and
#'   says so when it does.
#' @param rebuild_context Logical. Force a rebuild of the cached
#'   documentation context before asking, bypassing the normal
#'   version-matched cache. Rarely needed - see [RebuildSidekickContext()].
#'   Default `FALSE`.
#'
#' @return Invisibly, the assistant's reply (also printed).
#' @export
AskSidekick <- function(question = NULL, seurat_object = NULL, image_path = NULL,
                        same_dataset = NA, models = .sidekick_default_models,
                        rebuild_context = FALSE) {
  if (!requireNamespace("ellmer", quietly = TRUE))
    stop("Package 'ellmer' is required for AskSidekick(). Install with:\n",
         "  install.packages('ellmer')")

  # Capture the caller's literal expression for seurat_object (e.g. `Microglia`)
  # before it's touched, so generated code examples use the user's own
  # variable name instead of a generic placeholder like `seurat_obj`.
  obj_expr <- substitute(seurat_object)
  obj_var_name <- if (is.symbol(obj_expr)) deparse(obj_expr) else {
    txt <- paste(deparse(obj_expr), collapse = "")
    if (nzchar(txt) && nchar(txt) < 60) txt else "seurat_object"
  }

  if (is.null(question) && is.null(image_path)) {
    stop("Provide a question, an image_path, or both.")
  }
  key <- Sys.getenv("GEMINI_API_KEY")
  if (!nzchar(key))
    stop("GEMINI_API_KEY is not set. Get a free key at ",
         "https://aistudio.google.com/apikey, then add a line to ~/.Renviron:\n",
         "  GEMINI_API_KEY=your_key_here\nand restart R.")

  chunks <- .sidekick_get_context(rebuild = rebuild_context)

  extra_kw <- character(0)
  parts <- list()
  if (!is.null(image_path)) {
    if (!file.exists(image_path)) stop("Image not found: ", image_path)
    is_pdf <- grepl("\\.pdf$", image_path, ignore.case = TRUE)
    parts[[length(parts) + 1]] <- if (is_pdf)
      ellmer::content_pdf_file(image_path) else ellmer::content_image_file(image_path)

    # scSidekick's own output paths often ARE the naming convention - e.g.
    # RunSCssGSEA() creates a folder literally named "SC ssGSEA {label}", and
    # filenames like "...TrendLabeled.pdf" nearly spell out PlotTrendLabeled().
    # basename() alone misses the parent folder, which is often the stronger
    # signal. Use the last two path components, and feed their words into
    # retrieval too - a hardcoded keyword seed can't anticipate every
    # function name a path might mention.
    path_parts <- normalizePath(image_path, mustWork = FALSE)
    path_parts <- utils::tail(strsplit(path_parts, .Platform$file.sep)[[1]], 2)
    path_context <- paste(path_parts, collapse = "/")
    path_words <- unique(unlist(strsplit(path_context, "[^A-Za-z0-9]+")))
    path_words <- path_words[nchar(path_words) >= 3]

    extra_kw <- c(extra_kw, path_words,
                  "heatmap", "GroupHeatmap", "dotplot", "FastDotPlot",
                  "SplitDotPlot", "DimPlot", "PlotDimPlots", "spatial",
                  "PlotSpatialFeaturePlots", "volcano", "PlotVolcano",
                  "feature plot", "PlotFeaturePlots", "violin", "StackedVlnPlot",
                  "trend", "PlotTrendLabeled", "PlotTrendAndUMAP", "composition",
                  "GSEA", "ssGSEA", "pathway", "RunGSEA", "RunGSEA_pseudobulk",
                  "RunSCssGSEA", "RunEnrichment", "PlotGSEAEnrichment",
                  "PlotPathwayButterfly", "MetaSummary", "PlotMetaSummary",
                  "pseudotime", "ridge", "density", "trajectory",
                  "PlotPseudotime", "PlotFeatureDistribution", "PlotFeatureTrend",
                  "PlotTrajectory", "RunSlingshot", "correlation", "PlotCorrelation")

    # This guidance ALWAYS applies whenever an image is supplied - it must
    # not be skipped just because the user also typed a question (e.g. "I
    # want to recreate this figure").
    path_guidance <- paste0(
      "This image is a figure the user wants to recreate. Its path (folder + ",
      "filename) is: \"", path_context, "\".\n\n",
      "scSidekick functions hardcode their own output folder/file naming - it ",
      "is deterministic, not a style choice, and OVERRIDES any guess based on ",
      "how the figure looks. If the path contains a distinctive fragment that ",
      "matches a function's own naming convention (check the documentation ",
      "excerpts below for what each function actually names its outputs - ",
      "e.g. a folder starting \"SC ssGSEA\" is single-cell-level and made ONLY ",
      "by RunSCssGSEA(), never by RunGSEA_pseudobulk() even though both ",
      "compute ssGSEA scores; a filename containing \"TrendLabeled\" is made ",
      "ONLY by PlotTrendLabeled()), you MUST identify that exact function as ",
      "the source. Do not substitute a different function in the same family ",
      "because the visual layout looks more consistent with it - the path is ",
      "ground truth about which function actually ran; the visual layout is ",
      "not, since one function's output can look superficially similar to ",
      "another's.\n\n",
      "Using the path, the visible plot type/layout, and the user's ACTUAL ",
      "object metadata columns (given below): ",
      "(1) identify which scSidekick function most likely made this figure, ",
      "(2) for any grouping/coloring/faceting variable visible in the figure, ",
      "check whether a WORD IN THE FILENAME closely matches one of the user's ",
      "REAL column names (below) before falling back to a generic guess like ",
      "\"Cluster\" - prefer the real column that actually matches, ",
      "(3) READ any gene names, pathway names, or other row/column/axis text ",
      "labels actually visible in the image and transcribe them literally, ",
      "character for character - this is direct content you can see, not ",
      "something to infer or leave to the user to fill in later, ",
      "(4) give the R call to recreate something similar, using the labels ",
      "you transcribed in step 3 wherever the call needs specific gene or ",
      "pathway names rather than a generic placeholder. ",
      "When the user's own question below asks for something SPECIFIC visible ",
      "in the image (e.g. \"give me the genes shown\", \"make a heatmap with ",
      "these genes\", \"what pathways are in panel b\") rather than just ",
      "\"recreate this figure\", answer that literal request directly using ",
      "the transcribed labels from step 3 - do not substitute a generic ",
      "explanation of the pipeline that produced the whole figure for a ",
      "direct request to extract and reuse specific content from it. In that ",
      "case, prefer whichever scSidekick function most directly accepts that ",
      "content as a plain argument (e.g. a gene vector via `features = c(...)` ",
      "on a function like `GroupHeatmap()`) over the (possibly more complex) ",
      "function that produced the original figure, if that original function ",
      "instead requires prerequisite pipeline output the user hasn't said they ",
      "have (e.g. `gsea_files` pointing at already-saved GSEA result CSVs) - a ",
      "correct answer to \"build X from these genes\" does not route through a ",
      "function the user would first have to run something else to use. If you ",
      "build a gene/pathway vector from what you transcribed, that vector MUST ",
      "actually be passed into the function call as an argument - never define ",
      "it and then give a call that doesn't reference it. ",
      if (isTRUE(same_dataset)) paste0(
        "The user has CONFIRMED this figure was made from this exact object, ",
        "so there MUST be an exact matching column for every grouping/coloring ",
        "variable shown - search the full column list below carefully (objects ",
        "often have several similar-but-not-identical columns, e.g. `Cluster` ",
        "vs. `seurat_clusters` vs. `Assignment`) and use the one that exactly ",
        "matches, not just a plausible-sounding one. "
      ) else paste0(
        "The user has NOT confirmed this figure came from this exact object ",
        "(it may be from a different dataset, an earlier package version, or ",
        "a collaborator's file), so an exact column match may not exist. Give ",
        "your best estimate of which real column(s) are closest, and say ",
        "explicitly that it's an estimate rather than a confirmed match. "
      ),
      "If you're genuinely unsure which function or column applies, say so ",
      "explicitly and give your top 2-3 candidates with what would ",
      "distinguish them, rather than picking one with false confidence."
    )
    q <- if (is.null(question)) path_guidance
         else paste0(path_guidance, "\n\nThe user's own question: ", question)

    # Retrieval scoring gets a NARROW query (path + the user's own words
    # only) rather than the full path_guidance text above. path_guidance's
    # own instructional prose illustrates the naming-convention rule with
    # concrete example function names (RunSCssGSEA, PlotTrendLabeled, ...) -
    # those are meant to teach the MODEL the general pattern, but retrieval's
    # exact-function-mention bonus can't tell an illustrative example from a
    # real match, so the same fixed examples would win retrieval on every
    # single image-based question regardless of what the image actually is.
    retrieval_q <- paste(path_context, question %||% "")
  } else {
    q <- question
    retrieval_q <- q
  }

  q_for_object <- q
  if (!is.null(seurat_object)) {
    q_for_object <- paste0(q, "\n\n",
                           .sidekick_describe_object(seurat_object, obj_var_name))
  }

  system_prompt <- .sidekick_build_system_prompt(retrieval_q, chunks, extra_kw)

  parts[[length(parts) + 1]] <- q_for_object
  reply <- .sidekick_chat_with_fallback(system_prompt, parts, key, models)

  # Print explicitly rather than relying on echo from inside chat$chat()
  # (disabled above - see the comment there) or on R's own auto-print (the
  # return value is invisible(), by design, so a caller who does capture it
  # in a variable doesn't get it printed twice).
  cat(reply, "\n")

  gene_note <- .sidekick_check_gene_matches(reply, seurat_object, obj_var_name)
  if (!is.null(gene_note)) {
    reply <- paste0(reply, gene_note)
    message(gene_note)
  }

  col_note <- .sidekick_check_column_matches(reply, seurat_object, obj_var_name, same_dataset)
  if (!is.null(col_note)) {
    reply <- paste0(reply, col_note)
    message(col_note)
  }

  invisible(reply)
}
