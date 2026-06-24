# Internal plot helpers – shared between plot_beforeafter and plot_dist_beforeafter
# ponytail: deduplicates ~100 lines of theme/validation logic

#' Check a single positive numeric (or NULL)
#' @noRd
.check_positive_numeric <- function(val, name, allow_null = FALSE) {
  if (allow_null && is.null(val)) return(invisible(NULL))
  if (!is.numeric(val) || length(val) != 1L || is.na(val) || val <= 0)
    stop(sprintf("'%s' must be a single positive number.", name), call. = FALSE)
}

#' Resolve font size: explicit > global scale > base default
#' @noRd
.resolve_size <- function(explicit, scale_factor, base_default) {
  if (!is.null(explicit))     return(explicit)
  if (!is.null(scale_factor)) return(base_default * scale_factor)
  base_default
}

#' Build base ggplot2 theme from a theme name
#' @noRd
.build_base_theme <- function(theme_choice, base_size, font_family, strip_text_size = NULL) {
  valid_themes <- c("nature", "minimal", "classic", "bw", "light", "dark")
  theme_choice <- tolower(theme_choice)
  if (!theme_choice %in% valid_themes)
    stop(sprintf(
      "'theme' must be one of: %s. Got: '%s'.",
      paste(valid_themes, collapse = ", "), theme_choice
    ), call. = FALSE)

  base <- switch(
    theme_choice,
    "nature"  = ggplot2::theme_bw(base_size = base_size, base_family = font_family) +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_line(color = "grey92", linewidth = 0.4),
        panel.grid.minor = ggplot2::element_blank(),
        panel.border     = ggplot2::element_rect(color = "grey70", fill = NA, linewidth = 0.6),
        strip.background = ggplot2::element_rect(fill = "grey95", color = "grey70"),
        strip.text       = ggplot2::element_text(face = "bold", family = font_family,
                                                 size = strip_text_size)
      ),
    "minimal" = ggplot2::theme_minimal(base_size = base_size, base_family = font_family),
    "classic" = ggplot2::theme_classic(base_size = base_size, base_family = font_family),
    "bw"      = ggplot2::theme_bw(base_size     = base_size, base_family = font_family),
    "light"   = ggplot2::theme_light(base_size  = base_size, base_family = font_family),
    "dark"    = ggplot2::theme_dark(base_size   = base_size, base_family = font_family)
  )
  base
}

#' Build custom theme overlay with resolved font sizes
#' @noRd
.build_custom_theme <- function(font_family, title_size, xlab_size, ylab_size,
                                axis_text_size, legend_title_size = NULL,
                                legend_text_size = NULL, strip_text_size = NULL,
                                title_face = NULL) {
  th <- ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    plot.title       = ggplot2::element_text(size = title_size, hjust = 0.5,
                                             face = title_face, family = font_family),
    axis.title.x     = ggplot2::element_text(size = xlab_size,  family = font_family),
    axis.title.y     = ggplot2::element_text(size = ylab_size,  family = font_family),
    axis.text        = ggplot2::element_text(size = axis_text_size, family = font_family)
  )
  if (!is.null(legend_title_size))
    th <- th + ggplot2::theme(
      legend.title = ggplot2::element_text(size = legend_title_size, family = font_family)
    )
  if (!is.null(legend_text_size))
    th <- th + ggplot2::theme(
      legend.text = ggplot2::element_text(size = legend_text_size, family = font_family)
    )
  if (!is.null(strip_text_size))
    th <- th + ggplot2::theme(
      strip.text = ggplot2::element_text(size = strip_text_size, face = "bold",
                                         family = font_family)
    )
  th
}

#' Validate and resolve plot_what feature names
#' @noRd
.resolve_plot_what <- function(plot_what, available_cols, seed = NULL, max_default = 6L) {
  if (!is.null(plot_what)) {
    if (!is.character(plot_what) || length(plot_what) < 1L)
      stop("'plot_what' must be a non-empty character vector of feature names.", call. = FALSE)

    missing_features <- setdiff(plot_what, available_cols)
    if (length(missing_features) > 0L) {
      if (length(missing_features) == length(plot_what)) {
        warning(
          "Feature(s) not found in data: ",
          paste(missing_features, collapse = ", "),
          ". Falling back to random selection.",
          call. = FALSE
        )
        plot_what <- NULL
      } else {
        warning(
          "Feature(s) not found (skipped): ",
          paste(missing_features, collapse = ", "),
          call. = FALSE
        )
        plot_what <- intersect(plot_what, available_cols)
      }
    }
  }

  if (is.null(plot_what)) {
    if (length(available_cols) == 0L)
      stop("No features available to plot.", call. = FALSE)
    plot_what <- if (length(available_cols) <= max_default) {
      available_cols
    } else {
      sample(available_cols, max_default)
    }
    message(sprintf(
      "Automatically selected %d random features to plot: %s",
      length(plot_what), paste(plot_what, collapse = ", ")
    ))
  }
  plot_what
}

#' Save and restore RNG state around a block (avoids withr dependency)
#' @noRd
.with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  old_seed <- if (exists(".Random.seed", envir = globalenv()))
    get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    if (is.null(old_seed)) {
      rm(".Random.seed", envir = globalenv())
    } else {
      assign(".Random.seed", old_seed, envir = globalenv())
    }
  })
  set.seed(seed)
  expr
}
