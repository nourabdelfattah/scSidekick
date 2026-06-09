# =============================================================================
# scSidekick custom ggplot2 theme
# theme_NourMin: clean minimal theme with light panel borders and no gridlines.
# This is the canonical name for what appears as theme_min2() in analysis
# scripts - they are identical.
# =============================================================================

#' scSidekick minimal ggplot2 theme
#'
#' A clean, publication-ready theme built on [ggplot2::theme_light()].
#' Features: no gridlines, light grey panel borders, centered titles,
#' transparent strip backgrounds, and compact margins. Use this instead of
#' `theme_min2()` from analysis scripts - they are identical.
#'
#' @param base_size Base font size in points. Default `11`.
#' @param base_family Base font family. Default `""` (system default).
#' @return A ggplot2 theme object.
#' @export
theme_NourMin <- function(base_size = 11, base_family = "") {
  ggplot2::theme_light(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title        = ggplot2::element_text(size  = ggplot2::rel(0.9),
                                                hjust = 0.5, vjust = 0.5),
      panel.grid.major  = ggplot2::element_blank(),
      panel.grid.minor  = ggplot2::element_blank(),
      panel.background  = ggplot2::element_blank(),
      panel.border      = ggplot2::element_rect(fill   = NA,
                                                colour = "grey90",
                                                linewidth = 1),
      strip.background  = ggplot2::element_rect(fill = NA, colour = NA),
      strip.text.x      = ggplot2::element_text(colour = "black",
                                                size   = ggplot2::rel(1.2)),
      strip.text.y      = ggplot2::element_text(colour = "black",
                                                size   = ggplot2::rel(1.2)),
      title             = ggplot2::element_text(size  = ggplot2::rel(0.9),
                                                hjust = 0.5, vjust = 0.5),
      axis.text         = ggplot2::element_text(colour = "black",
                                                size   = ggplot2::rel(0.8)),
      axis.title        = ggplot2::element_blank(),
      legend.title      = ggplot2::element_text(colour = "black",
                                                size   = ggplot2::rel(0.9),
                                                hjust  = 0.5),
      legend.key.size   = ggplot2::unit(0.9, "lines"),
      legend.text       = ggplot2::element_text(size   = ggplot2::rel(0.7),
                                                colour = "black"),
      legend.key        = ggplot2::element_rect(colour = NA, fill = NA),
      legend.background = ggplot2::element_rect(colour = NA, fill = NA),
      plot.margin       = ggplot2::unit(c(0.1, 0, 0, -0.2), "lines")
    )
}
