# =============================================================================
# scSidekick CellChat visualisation helpers
#
# These are lightly adapted versions of the "NetViasualHack" functions
# originally written for the Yun lab by Nhat Nguyen (Box-Box repository).
# The key difference from standard CellChat visualisation functions is that
# they access cc@net$prob (individual L-R pair probabilities, pre-aggregation)
# rather than cc@netP$prob (pathway-level significance-filtered aggregates).
# A user-controlled `thresh` parameter (p-value cutoff on cc@net$pval) lets
# callers decide how strict to be: thresh = 0.05 matches CellChat's default
# behaviour; thresh = 1 shows ALL predicted interactions regardless of
# significance, enabling comparison of pathways that exist in some conditions
# but did not reach pathway-level significance in others.
#
# Functions:
#   .cc_search_pair          - inline replacement for CellChat:::searchPair
#   .cc_chord_internal       - chord diagram with group-merge support
#   .cc_aggregate            - aggregate visualisation (circle / chord / hierarchy)
#   .cc_chord_cell           - pathway chord with slot.name choice
#   .cc_heatmap              - heatmap with zero-matrix fallback for absent paths
# =============================================================================

# ---------------------------------------------------------------------------
# .cc_search_pair  - find L-R rows for a given pathway in an LRsig table
# ---------------------------------------------------------------------------
.cc_search_pair <- function(signaling, lrsig) {
  lrsig[!is.na(lrsig$pathway_name) &
        lrsig$pathway_name == signaling, , drop = FALSE]
}

# ---------------------------------------------------------------------------
# .cc_chord_internal
# Modified CellChat chord-cell function that forwards the `group` parameter
# to circlize::chordDiagram for cell-type merging.
# ---------------------------------------------------------------------------
.cc_chord_internal <- function(
    net,
    color.use            = NULL,
    group                = NULL,
    cell.order           = NULL,
    sources.use          = NULL,
    targets.use          = NULL,
    lab.cex              = 0.8,
    small.gap            = 1,
    big.gap              = 10,
    annotationTrackHeight = c(0.03),
    remove.isolate       = FALSE,
    link.visible         = TRUE,
    scale                = FALSE,
    directional          = 1,
    link.target.prop     = TRUE,
    reduce               = -1,
    transparency         = 0.4,
    link.border          = NA,
    title.name           = NULL,
    show.legend          = FALSE,
    legend.pos.x         = 20,
    legend.pos.y         = 20,
    ...
) {
  if (inherits(net, c("matrix", "Matrix"))) {
    cell.levels <- union(rownames(net), colnames(net))
    net         <- reshape2::melt(net, value.name = "prob")
    colnames(net)[1:2] <- c("source", "target")
  } else if (is.data.frame(net)) {
    if (!all(c("source", "target", "prob") %in% colnames(net)))
      stop("net must have columns: source, target, prob")
    cell.levels <- as.character(union(net$source, net$target))
  }

  if (!is.null(cell.order)) cell.levels <- cell.order
  net$source <- as.character(net$source)
  net$target <- as.character(net$target)

  if (!is.null(sources.use)) {
    if (is.numeric(sources.use)) sources.use <- cell.levels[sources.use]
    net <- subset(net, source %in% sources.use)
  }
  if (!is.null(targets.use)) {
    if (is.numeric(targets.use)) targets.use <- cell.levels[targets.use]
    net <- subset(net, target %in% targets.use)
  }

  net <- subset(net, prob > 0)
  if (nrow(net) == 0) message("No interactions to display for this pathway.")

  if (!remove.isolate) {
    cells.removed <- setdiff(cell.levels,
                             as.character(union(net$source, net$target)))
    if (length(cells.removed) > 0) {
      net.fake <- data.frame(
        cells.removed, cells.removed,
        1e-10 * sample(length(cells.removed), length(cells.removed))
      )
      colnames(net.fake) <- colnames(net)
      net          <- rbind(net, net.fake)
      link.visible <- net[, 1:2]
      link.visible$plot <- FALSE
      if (nrow(net) > nrow(net.fake))
        link.visible$plot[seq_len(nrow(net) - nrow(net.fake))] <- TRUE
      scale <- TRUE
    }
  }

  cells.use    <- union(net$source, net$target)
  order.sector <- cell.levels[cell.levels %in% cells.use]

  if (is.null(color.use)) {
    color.use          <- Nour_pal("all")(length(cell.levels))
    names(color.use)   <- cell.levels
  } else if (is.null(names(color.use))) {
    names(color.use)   <- cell.levels
  }

  grid.col        <- color.use[order.sector]
  names(grid.col) <- order.sector

  if (!is.null(group))
    group <- group[names(group) %in% order.sector]

  edge.color  <- color.use[as.character(net$source)]
  link.arr.type <- if (directional %in% c(0, 2)) "triangle" else "big.arrow"

  circlize::circos.clear()
  circlize::chordDiagram(
    net,
    order               = order.sector,
    col                 = edge.color,
    grid.col            = grid.col,
    transparency        = transparency,
    link.border         = link.border,
    directional         = directional,
    direction.type      = c("diffHeight", "arrows"),
    link.arr.type       = link.arr.type,
    annotationTrack     = "grid",
    annotationTrackHeight = annotationTrackHeight,
    preAllocateTracks   = list(track.height = max(strwidth(order.sector))),
    small.gap           = small.gap,
    big.gap             = big.gap,
    link.visible        = link.visible,
    scale               = scale,
    group               = group,
    link.target.prop    = link.target.prop,
    reduce              = reduce,
    ...
  )

  circlize::circos.track(track.index = 1, panel.fun = function(x, y) {
    xlim        <- circlize::get.cell.meta.data("xlim")
    ylim        <- circlize::get.cell.meta.data("ylim")
    sector.name <- circlize::get.cell.meta.data("sector.index")
    circlize::circos.text(mean(xlim), ylim[1], sector.name,
                          facing = "clockwise", niceFacing = TRUE,
                          adj = c(0, 0.5), cex = lab.cex)
  }, bg.border = NA)

  if (show.legend) {
    lgd <- ComplexHeatmap::Legend(
      at        = names(grid.col),
      type      = "grid",
      legend_gp = grid::gpar(fill = grid.col),
      title     = "Cell State"
    )
    ComplexHeatmap::draw(
      lgd,
      x    = grid::unit(1, "npc") - grid::unit(legend.pos.x, "mm"),
      y    = grid::unit(legend.pos.y, "mm"),
      just = c("right", "bottom")
    )
  }
  if (!is.null(title.name))
    graphics::text(-0, 1.02, title.name, cex = 2)

  circlize::circos.clear()
  invisible(grDevices::recordPlot())
}

# ---------------------------------------------------------------------------
# .cc_aggregate
# Like CellChat::netVisual_aggregate but reads cc@net$prob so it can display
# pathways that have L-R interactions even when not significant at the
# pathway level.  thresh controls which L-R pairs are shown (0.05 = default
# CellChat significance; 1 = show everything).
# ---------------------------------------------------------------------------
.cc_aggregate <- function(
    object,
    signaling,
    signaling.name  = NULL,
    color.use       = NULL,
    thresh          = 0.05,
    vertex.receiver = NULL,
    sources.use     = NULL,
    targets.use     = NULL,
    idents.use      = NULL,
    remove.isolate  = FALSE,
    vertex.weight   = NULL,
    vertex.size.max = NULL,
    weight.scale    = TRUE,
    edge.weight.max = NULL,
    edge.width.max  = 8,
    layout          = c("circle", "chord", "hierarchy"),
    group           = NULL,
    cell.order      = NULL,
    small.gap       = 1,
    big.gap         = 10,
    scale           = FALSE,
    reduce          = -1,
    show.legend     = FALSE,
    legend.pos.x    = 20,
    legend.pos.y    = 20,
    ...
) {
  layout <- match.arg(layout)

  if (is.null(vertex.weight))
    vertex.weight <- as.numeric(table(object@idents))
  if (is.null(vertex.size.max))
    vertex.size.max <- if (length(unique(vertex.weight)) == 1) 5 else 15
  if (is.null(signaling.name))
    signaling.name <- signaling

  # Find L-R pairs for this pathway in the object's LRsig
  pairLR      <- .cc_search_pair(signaling, object@LR$LRsig)
  pairLR.names <- intersect(rownames(pairLR), dimnames(object@net$prob)[[3]])
  if (length(pairLR.names) == 0)
    stop("No L-R pairs found for pathway '", signaling, "'.")

  pairLR <- pairLR[pairLR.names, , drop = FALSE]
  prob   <- object@net$prob[, , pairLR.names, drop = FALSE]
  pval   <- object@net$pval[, , pairLR.names, drop = FALSE]

  # Apply threshold - zero out non-significant interactions
  prob[pval > thresh] <- 0

  # Ensure 3-D even for a single L-R pair
  if (length(dim(prob)) == 2)
    prob <- array(prob, dim = c(dim(prob), 1),
                  dimnames = c(dimnames(prob), list(pairLR.names)))

  prob.sum <- apply(prob, c(1, 2), sum)

  if (layout == "circle") {
    if (is.null(edge.weight.max)) edge.weight.max <- max(prob.sum)
    CellChat::netVisual_circle(
      prob.sum,
      sources.use     = sources.use,
      targets.use     = targets.use,
      idents.use      = idents.use,
      remove.isolate  = remove.isolate,
      color.use       = color.use,
      vertex.weight   = vertex.weight,
      vertex.size.max = vertex.size.max,
      weight.scale    = weight.scale,
      edge.weight.max = edge.weight.max,
      edge.width.max  = edge.width.max,
      title.name      = paste0(signaling.name, " signaling pathway network"),
      ...
    )
  } else if (layout == "chord") {
    .cc_chord_internal(
      prob.sum,
      color.use    = color.use,
      sources.use  = sources.use,
      targets.use  = targets.use,
      remove.isolate = remove.isolate,
      group        = group,
      cell.order   = cell.order,
      small.gap    = small.gap,
      big.gap      = big.gap,
      scale        = scale,
      reduce       = reduce,
      title.name   = paste0(signaling.name, " signaling pathway network"),
      show.legend  = show.legend,
      legend.pos.x = legend.pos.x,
      legend.pos.y = legend.pos.y
    )
  } else if (layout == "hierarchy") {
    if (is.null(edge.weight.max)) edge.weight.max <- max(prob.sum)
    graphics::par(mfrow = c(1, 2))
    CellChat::netVisual_hierarchy1(
      prob.sum,
      vertex.receiver = vertex.receiver,
      sources.use     = sources.use,
      targets.use     = targets.use,
      remove.isolate  = remove.isolate,
      color.use       = color.use,
      vertex.weight   = vertex.weight,
      vertex.size.max = vertex.size.max,
      weight.scale    = weight.scale,
      edge.weight.max = edge.weight.max,
      edge.width.max  = edge.width.max,
      ...
    )
    CellChat::netVisual_hierarchy2(
      prob.sum,
      vertex.receiver = setdiff(seq_len(nrow(prob.sum)), vertex.receiver),
      sources.use     = sources.use,
      targets.use     = targets.use,
      remove.isolate  = remove.isolate,
      color.use       = color.use,
      vertex.weight   = vertex.weight,
      vertex.size.max = vertex.size.max,
      weight.scale    = weight.scale,
      edge.weight.max = edge.weight.max,
      edge.width.max  = edge.width.max,
      ...
    )
    graphics::mtext(paste0(signaling.name, " signaling pathway network"),
                    side = 3, outer = TRUE, cex = 1, line = -6)
  }

  invisible(grDevices::recordPlot())
}

# ---------------------------------------------------------------------------
# .cc_chord_cell
# Like netVisual_chord_cell2: pathway chord diagram with group-merge support
# and threshold-controlled significance filtering.
# ---------------------------------------------------------------------------
.cc_chord_cell <- function(
    object,
    signaling     = NULL,
    slot.name     = "netP",
    color.use     = NULL,
    group         = NULL,
    cell.order    = NULL,
    sources.use   = NULL,
    targets.use   = NULL,
    lab.cex       = 0.8,
    small.gap     = 1,
    big.gap       = 10,
    annotationTrackHeight = c(0.03),
    remove.isolate = FALSE,
    link.visible   = TRUE,
    scale          = FALSE,
    directional    = 1,
    link.target.prop = TRUE,
    reduce         = -1,
    transparency   = 0.4,
    link.border    = NA,
    title.name     = NULL,
    show.legend    = FALSE,
    legend.pos.x   = 20,
    legend.pos.y   = 20,
    thresh         = 0.05,
    ...
) {
  if (is.null(signaling)) stop("Please provide `signaling`.")

  pairLR       <- .cc_search_pair(signaling, object@LR$LRsig)
  pairLR.names <- intersect(rownames(pairLR), dimnames(object@net$prob)[[3]])
  if (length(pairLR.names) == 0)
    stop("No L-R pairs found for pathway '", signaling, "'.")

  prob <- object@net$prob[, , pairLR.names, drop = FALSE]
  pval <- object@net$pval[, , pairLR.names, drop = FALSE]
  prob[pval > thresh] <- 0

  if (length(dim(prob)) == 2)
    prob <- array(prob, dim = c(dim(prob), 1))

  net <- apply(prob, c(1, 2), sum)
  if (is.null(title.name))
    title.name <- paste0(signaling, " signaling pathway network")

  .cc_chord_internal(
    net,
    color.use    = color.use,
    group        = group,
    cell.order   = cell.order,
    sources.use  = sources.use,
    targets.use  = targets.use,
    lab.cex      = lab.cex,
    small.gap    = small.gap,
    big.gap      = big.gap,
    annotationTrackHeight = annotationTrackHeight,
    remove.isolate = remove.isolate,
    link.visible = link.visible,
    scale        = scale,
    directional  = directional,
    link.target.prop = link.target.prop,
    reduce       = reduce,
    transparency = transparency,
    link.border  = link.border,
    title.name   = title.name,
    show.legend  = show.legend,
    legend.pos.x = legend.pos.x,
    legend.pos.y = legend.pos.y,
    ...
  )
}

# ---------------------------------------------------------------------------
# .cc_heatmap
# Like netVisual_heatmap2: reads netP$prob[,,signaling] with a zero-matrix
# fallback when the pathway is absent from netP (i.e., not significant at
# the pathway level but present at the L-R level).
# ---------------------------------------------------------------------------
# Extra arguments in ... are forwarded to ComplexHeatmap::Heatmap(), so users
# can override any default - e.g. column_names_rot, row_dend_width, col, etc.
.cc_heatmap <- function(
    object,
    signaling       = NULL,
    color.use       = NULL,
    color.heatmap   = "Reds",
    title.name      = NULL,
    font.size       = 8,
    font.size.title = 10,
    cluster.rows    = FALSE,
    cluster.cols    = FALSE,
    sources.use     = NULL,
    targets.use     = NULL,
    remove.isolate  = FALSE,
    ...
) {
  # Try pathway-level probability (netP); fall back to zero matrix
  net.diff <- tryCatch(
    object@netP$prob[, , signaling, drop = TRUE],
    error = function(e) {
      # Pathway not in netP - create zero matrix of same dims as net$prob
      n <- nrow(object@net$prob)
      m <- matrix(0, n, n,
                  dimnames = list(rownames(object@net$prob),
                                  colnames(object@net$prob)))
      m
    }
  )

  if (is.null(title.name))
    title.name <- paste0(signaling, " signaling network")

  legend.name <- "Communication Prob."

  if (!is.null(sources.use) || !is.null(targets.use)) {
    df.net          <- reshape2::melt(net.diff, value.name = "value")
    colnames(df.net)[1:2] <- c("source", "target")
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) sources.use <- rownames(net.diff)[sources.use]
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) targets.use <- rownames(net.diff)[targets.use]
      df.net <- subset(df.net, target %in% targets.use)
    }
    df.net$source <- factor(df.net$source, levels = rownames(net.diff))
    df.net$target <- factor(df.net$target, levels = rownames(net.diff))
    df.net$value[is.na(df.net$value)] <- 0
    net.diff <- tapply(df.net$value, list(df.net$source, df.net$target), sum)
  }
  net.diff[is.na(net.diff)] <- 0

  if (remove.isolate) {
    idx <- intersect(which(rowSums(net.diff) == 0),
                     which(colSums(net.diff) == 0))
    if (length(idx) > 0) { net.diff <- net.diff[-idx, -idx] }
  }

  mat <- net.diff

  if (is.null(color.use)) {
    color.use <- Nour_pal("all")(ncol(mat))
  }
  names(color.use) <- colnames(mat)

  # Colour function
  if (sum(abs(mat)) == 0) {
    col_fun <- "white"
  } else if (length(color.heatmap) == 1) {
    col_fun <- grDevices::colorRampPalette(
      RColorBrewer::brewer.pal(9, color.heatmap))(100)
  } else {
    col_fun <- circlize::colorRamp2(
      c(min(mat), max(mat)), color.heatmap)
  }

  df_ann  <- data.frame(group = colnames(mat), row.names = colnames(mat))
  col_ann <- ComplexHeatmap::HeatmapAnnotation(
    df = df_ann,
    col = list(group = color.use),
    which = "column",
    show_legend = FALSE, show_annotation_name = FALSE,
    simple_anno_size = grid::unit(0.2, "cm")
  )
  row_ann <- ComplexHeatmap::HeatmapAnnotation(
    df = df_ann,
    col = list(group = color.use),
    which = "row",
    show_legend = FALSE, show_annotation_name = FALSE,
    simple_anno_size = grid::unit(0.2, "cm")
  )
  ha1 <- ComplexHeatmap::rowAnnotation(
    Strength = ComplexHeatmap::anno_barplot(
      rowSums(abs(mat)), border = FALSE,
      gp = grid::gpar(fill = color.use, col = color.use)),
    show_annotation_name = FALSE
  )
  ha2 <- ComplexHeatmap::HeatmapAnnotation(
    Strength = ComplexHeatmap::anno_barplot(
      colSums(abs(mat)), border = FALSE,
      gp = grid::gpar(fill = color.use, col = color.use)),
    show_annotation_name = FALSE
  )

  if (sum(abs(mat)) == 0)
    mat_plot <- mat
  else {
    mat_plot <- mat
    mat_plot[mat_plot == 0] <- NA
  }

  # Build argument list; merge with ... so callers can override any default
  ht_args <- utils::modifyList(
    list(
      mat_plot,
      col                  = col_fun,
      na_col               = "white",
      name                 = legend.name,
      bottom_annotation    = col_ann,
      left_annotation      = row_ann,
      top_annotation       = ha2,
      right_annotation     = ha1,
      cluster_rows         = cluster.rows,
      cluster_columns      = cluster.cols,
      row_names_side       = "left",
      row_names_rot        = 0,
      row_names_gp         = grid::gpar(fontsize = font.size),
      column_names_gp      = grid::gpar(fontsize = font.size),
      column_title         = title.name,
      column_title_gp      = grid::gpar(fontsize = font.size.title),
      column_names_rot     = 90,
      row_title            = "Sources (Sender)",
      row_title_gp         = grid::gpar(fontsize = font.size.title),
      row_title_rot        = 90,
      heatmap_legend_param = list(
        title_gp        = grid::gpar(fontsize = 8, fontface = "plain"),
        title_position  = "leftcenter-rot",
        border          = NA,
        legend_height   = grid::unit(20, "mm"),
        labels_gp       = grid::gpar(fontsize = 8),
        grid_width      = grid::unit(2, "mm")
      )
    ),
    list(...)
  )
  do.call(ComplexHeatmap::Heatmap, ht_args)
}
