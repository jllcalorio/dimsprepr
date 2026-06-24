#' Plot Pairwise Before-and-After Comparison for Selected Features
#'
#' @title Plot Pairwise Before-and-After Comparison for Selected Features
#'
#' @description
#' Generates a multi-panel plot comparing the signal of selected features
#' before and after any data transformation, correction, or manipulation.
#' Each panel displays one feature across the injection sequence, colored by
#' sample group and shaped by batch, with an optional LOESS trend line fitted
#' to a user-defined reference group (e.g., QC samples).
#'
#' @param x A \code{data.frame}, \code{tibble}, or \code{matrix} representing
#'   the data \emph{before} the transformation. Rows are samples; columns are
#'   features. Column names may contain special characters. Input can be
#'   a \code{link{run_DIpreprocess}} class.
#' @param y A \code{data.frame}, \code{tibble}, or \code{matrix} representing
#'   the data \emph{after} the transformation. Must have the same dimensions
#'   and column names as \code{x}.
#' @param metadata A \code{data.frame} containing sample-level metadata shared
#'   by both \code{x} and \code{y}. Must have the same number of rows as
#'   \code{x} and \code{y}.
#' @param plot_what Character vector or \code{NULL}. Column names from \code{x} to be plotted.
#'   If \code{NULL} (default), the function randomly selects up to 6 features.
#'   If specific features are provided but not found, a warning is issued and they
#'   are skipped. If none of the provided features are found, the function falls
#'   back to random selection.
#' @param col_injection Character. Name of the column in \code{metadata}
#'   containing the injection sequence (numeric order of runs).
#'   Default is \code{"InjectionSequence"}. If \code{run_DIpreprocess} class,
#'   value is automatically taken from `$parameters$injection_col`.
#' @param col_batch Character. Name of the column in \code{metadata}
#'   containing batch identifiers. Default is \code{"Batch"}. If \code{run_DIpreprocess} class,
#'   value is automatically taken from `$parameters$batch_col`.
#' @param col_group Character. Name of the column in \code{metadata}
#'   containing sample group labels (e.g., \code{"QC"}, \code{"Sample"}).
#'   Default is \code{"Group"}. If \code{run_DIpreprocess} class,
#'   value is automatically taken from `$parameters$group_col`.
#' @param col_qc_label Character vector. The value(s) in the \code{col_group} column
#'   that identifies the reference group for LOESS trend fitting (typically
#'   QC samples). Can be a vector like \code{c("SQC", "EQC")}. If \code{NULL},
#'   no LOESS line is drawn. Default is \code{c("QC", "SQC", "EQC")}. If \code{run_DIpreprocess} class,
#'   value is automatically taken from `$parameters$qc_types`.
#' @param arrange_group Character vector or \code{NULL}. Custom order for the
#'   sample groups in \code{col_group}. If \code{NULL} (default), groups are
#'   arranged alphanumerically, with \code{col_qc_label} values placed first.
#' @param plot_title Character or \code{NULL}. Main title for the combined plot.
#'   If \code{NULL} (default), a title is automatically generated based on
#'   \code{x_label}, \code{y_label}, and the number of batches.
#' @param x_label Character. Label for the before-correction facet panel.
#'   Default is \code{"Before"}.
#' @param y_label Character. Label for the after-correction facet panel.
#'   Default is \code{"After"}.
#' @param log_transform Logical. If \code{TRUE}, applies a \code{log10}
#'   transformation (with a small floor of \code{1e-10}) to feature values
#'   before plotting. Default is \code{TRUE}.
#' @param point_size Numeric. Size of individual data points. Default is
#'   \code{2}.
#' @param point_alpha Numeric in \code{[0, 1]}. Transparency of data points.
#'   Default is \code{0.7}.
#' @param plot_cols Integer or \code{NULL}. Number of columns in the combined
#'   plot grid. If \code{NULL}, determined automatically. Default is
#'   \code{NULL}.
#' @param theme Character. ggplot2 theme to apply to each panel. Options:
#'   \code{"nature"}, \code{"minimal"}, \code{"classic"}, \code{"bw"},
#'   \code{"light"}, \code{"dark"}. Default is \code{"nature"} (a clean,
#'   publication-ready style based on \code{theme_bw()}).
#' @param base_size Numeric. Base font size for the theme (pts). Default is
#'   \code{15}.
#' @param font_family Character. Font family for all text elements. Default
#'   is \code{"sans"}.
#' @param seed Numeric or \code{NULL}. Random seed passed to \code{set.seed()}
#'   before any internal sampling operations. Default is \code{123}.
#' @param global_font_size Numeric or \code{NULL}. Base font size for
#'   proportional scaling of all text elements. Default is \code{NULL}.
#' @param title_size Numeric or \code{NULL}. Feature panel title font size.
#'   Default is \code{NULL}.
#' @param subtitle_size Numeric or \code{NULL}. Main plot title (top grob)
#'   font size. Default is \code{NULL}.
#' @param xlab_size Numeric or \code{NULL}. X-axis label font size. Default
#'   is \code{NULL}.
#' @param ylab_size Numeric or \code{NULL}. Y-axis label font size. Default
#'   is \code{NULL}.
#'
#' @details
#' \code{x} and \code{y} must share identical row counts, column counts, and
#' column names. Row order is assumed to correspond between \code{x},
#' \code{y}, and \code{metadata}.
#'
#' When \code{log_transform = TRUE}, values are floored at \code{1e-10}
#' before log transformation to avoid \code{-Inf} results from zero or
#' negative values.
#'
#' Font sizes are resolved with the following priority (highest to lowest):
#' explicit size argument (e.g., \code{title_size}) > \code{global_font_size}
#' scaling > \code{base_size}-derived default. When \code{global_font_size}
#' is provided and a specific size is not, text elements scale proportionally
#' relative to a base size of 11.
#'
#' The LOESS trend line is fitted only to the group identified by
#' \code{col_qc_label}. If fewer than 3 observations of that group are
#' present, a warning is issued and the line is suppressed for that feature.
#'
#' The combined plot is rendered via \code{gridExtra::arrangeGrob()} and
#' displayed with \code{grid::grid.draw()}. The returned object is an
#' invisible \code{gtable} grob that can be redrawn or exported with
#' \code{ggplot2::ggsave()}.
#'
#' @return Invisibly returns a \code{gtable} grob object (from
#'   \code{gridExtra::arrangeGrob()}) containing all feature panels. The
#'   plot is also drawn to the active graphics device.
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' n_samples  <- 30
#' n_features <- 4
#'
#' x <- as.data.frame(
#'   matrix(abs(rnorm(n_samples * n_features, mean = 1000, sd = 200)),
#'          nrow = n_samples,
#'          dimnames = list(NULL, paste0("Feature_", seq_len(n_features))))
#' )
#' y <- as.data.frame(x + rnorm(n_samples * n_features, sd = 50))
#'
#' meta <- data.frame(
#'   InjectionSequence = seq_len(n_samples),
#'   Batch             = rep(c("B1", "B2"), each = n_samples / 2),
#'   Group             = rep(c("Sample", "QC"), times = c(n_samples - 6, 6))
#' )
#'
#' plot_beforeafter(
#'   x             = x,
#'   y             = y,
#'   metadata      = meta,
#'   plot_what     = c("Feature_1", "Feature_2"),
#'   col_injection = "InjectionSequence",
#'   col_batch     = "Batch",
#'   col_group     = "Group",
#'   col_qc_label  = "QC",
#'   theme         = "nature"
#' )
#' }
#'
#' @author John Lennon L. Calorio
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth facet_wrap
#' @importFrom ggplot2 scale_color_manual scale_shape_manual labs theme_bw theme_minimal theme_classic theme_light
#' @importFrom ggplot2 theme_dark theme element_text element_blank element_line element_rect
#' @importFrom ggplot2 labeller
#' @importFrom gridExtra arrangeGrob
#' @importFrom grid grid.newpage grid.draw textGrob gpar
#'
#' @export
plot_beforeafter <- function(
    x,
    y,
    metadata,
    plot_what        = NULL,
    col_injection    = "InjectionSequence",
    col_batch        = "Batch",
    col_group        = "Group",
    col_qc_label     = c("QC", "SQC", "EQC"),
    arrange_group    = NULL,
    plot_title       = NULL,
    x_label          = "Before",
    y_label          = "After",
    log_transform    = TRUE,
    point_size       = 2,
    point_alpha      = 0.7,
    plot_cols        = NULL,
    theme            = "nature",
    base_size        = 15,
    font_family      = "sans",
    seed             = 123,
    global_font_size = NULL,
    title_size       = NULL,
    subtitle_size    = NULL,
    xlab_size        = NULL,
    ylab_size        = NULL
) {

  # --- Integration with run_DIpreprocess ---
  if (inherits(x, "run_DIpreprocess")) {
    obj_ref <- x
    if (missing(metadata) || is.null(metadata)) {
      # Must use unmerged metadata because data_imputed is unmerged
      metadata <- x$metadata
    }
    
    # Extract states
    data_before <- x$data_imputed
    data_after  <- x$data_nonpls
    
    # Subset to common features to ensure exact alignment
    common_features <- intersect(colnames(data_before), colnames(data_after))
    if (length(common_features) == 0) {
      stop("No common features found between imputed and final non-PLS data.")
    }
    
    # Overwrite x and y
    x <- data_before[, common_features, drop = FALSE]
    y <- data_after[, common_features, drop = FALSE]

    # Use standardized column names from run_DIpreprocess metadata
    col_injection <- "InjectionSequence"
    col_batch     <- "Batch"
    col_group     <- "Group"

    # Use qc_types from the parameters for LOESS reference
    if (!is.null(obj_ref$parameters$qc_types)) {
      # Evaluate because parameters might be stored as calls/symbols
      col_qc_label <- eval(obj_ref$parameters$qc_types, envir = parent.frame())
    } else {
      col_qc_label <- c("QC", "SQC", "EQC") # Fallback to defaults
    }
  }

  # -------------------------------------------------------------------------
  # Input validation: x and y
  # -------------------------------------------------------------------------
  .is_tabular <- function(z) is.data.frame(z) || is.matrix(z)

  if (!.is_tabular(x)) stop("'x' must be a data.frame, tibble, or matrix.")
  if (!.is_tabular(y)) stop("'y' must be a data.frame, tibble, or matrix.")

  if (!identical(dim(x), dim(y))) {
    stop(
      "'x' and 'y' must have the same dimensions. ",
      "Got: x = ", nrow(x), " x ", ncol(x), ", ",
      "y = ",      nrow(y), " x ", ncol(y), "."
    )
  }

  x_cols <- if (is.matrix(x)) colnames(x) else names(x)
  y_cols <- if (is.matrix(y)) colnames(y) else names(y)

  if (!identical(x_cols, y_cols)) {
    stop("'x' and 'y' must have identical column names in the same order.")
  }
  if (is.null(x_cols)) stop("'x' and 'y' must have named columns.")

  x <- as.data.frame(x, check.names = FALSE)
  y <- as.data.frame(y, check.names = FALSE)

  # -------------------------------------------------------------------------
  # Input validation: metadata
  # -------------------------------------------------------------------------
  if (!is.data.frame(metadata)) stop("'metadata' must be a data.frame.")
  if (nrow(metadata) != nrow(x)) {
    stop(
      "'metadata' must have the same number of rows as 'x' and 'y'. ",
      "Got: metadata = ", nrow(metadata), " rows, x = ", nrow(x), " rows."
    )
  }

  for (col_arg in list(
    list(val = col_injection, name = "col_injection"),
    list(val = col_batch,     name = "col_batch"),
    list(val = col_group,     name = "col_group")
  )) {
    if (!is.character(col_arg$val) || length(col_arg$val) != 1L) {
      stop("'", col_arg$name, "' must be a single character string.")
    }
    if (!col_arg$val %in% names(metadata)) {
      stop(
        "Column '", col_arg$val, "' not found in 'metadata'. ",
        "Available columns: ", paste(names(metadata), collapse = ", "), "."
      )
    }
  }

  # -------------------------------------------------------------------------
  # Validate and Resolve plot_what
  # -------------------------------------------------------------------------
  if (!is.null(plot_what)) {
    if (!is.character(plot_what) || length(plot_what) < 1L) {
      stop("'plot_what' must be a non-empty character vector of feature names.")
    }

    missing_features <- setdiff(plot_what, x_cols)
    if (length(missing_features) > 0L) {
      if (length(missing_features) == length(plot_what)) {
        warning(
          "The following feature(s) in 'plot_what' were not found in 'x'/'y': ",
          paste(missing_features, collapse = ", "), ". ",
          "Falling back to random selection.",
          call. = FALSE
        )
        plot_what <- NULL
      } else {
        warning(
          "The following feature(s) in 'plot_what' were not found in 'x'/'y' and will be skipped: ",
          paste(missing_features, collapse = ", "), ".",
          call. = FALSE
        )
        plot_what <- intersect(plot_what, x_cols)
      }
    }
  }

  if (is.null(plot_what)) {
    if (length(x_cols) == 0L) stop("No features available in 'x' to plot.")

    plot_what <- .with_seed(seed, {
      if (length(x_cols) <= 6) x_cols else sample(x_cols, 6)
    })

    message(sprintf(
      "Automatically selected %d random features to plot: %s",
      length(plot_what), paste(plot_what, collapse = ", ")
    ))
  }

  # -------------------------------------------------------------------------
  # Input validation: scalar numerics
  # -------------------------------------------------------------------------
  .check_positive_numeric(point_size,       "point_size")
  .check_positive_numeric(point_alpha,      "point_alpha")
  .check_positive_numeric(base_size,        "base_size")
  .check_positive_numeric(global_font_size, "global_font_size", allow_null = TRUE)
  .check_positive_numeric(title_size,       "title_size",       allow_null = TRUE)
  .check_positive_numeric(subtitle_size,    "subtitle_size",    allow_null = TRUE)
  .check_positive_numeric(xlab_size,        "xlab_size",        allow_null = TRUE)
  .check_positive_numeric(ylab_size,        "ylab_size",        allow_null = TRUE)

  if (point_alpha > 1) {
    warning("'point_alpha' is > 1; it will be clamped to 1.")
    point_alpha <- 1
  }

  if (!is.null(plot_cols)) {
    if (!is.numeric(plot_cols) || length(plot_cols) != 1L ||
        is.na(plot_cols) || plot_cols < 1L || plot_cols != floor(plot_cols)) {
      stop("'plot_cols' must be a single positive integer or NULL.")
    }
    plot_cols <- as.integer(plot_cols)
  }

  if (!is.logical(log_transform) || length(log_transform) != 1L) {
    stop("'log_transform' must be a single logical value (TRUE or FALSE).")
  }

  # -------------------------------------------------------------------------
  # Input validation: theme and font
  # -------------------------------------------------------------------------
  valid_themes <- c("nature", "minimal", "classic", "bw", "light", "dark")
  theme_choice <- tolower(theme)
  if (!theme_choice %in% valid_themes) {
    stop(
      "'theme' must be one of: ", paste(valid_themes, collapse = ", "), ". ",
      "Got: '", theme, "'."
    )
  }

  if (!is.character(font_family) || length(font_family) != 1L) {
    stop("'font_family' must be a single character string.")
  }

  # -------------------------------------------------------------------------
  # Input validation: col_qc_label
  # -------------------------------------------------------------------------
  if (!is.null(col_qc_label)) {
    if (!is.character(col_qc_label) || length(col_qc_label) < 1L) {
      stop("'col_qc_label' must be a character vector or NULL.")
    }
    if (!any(col_qc_label %in% metadata[[col_group]])) {
      warning(
        "None of the 'col_qc_label' values (", paste(col_qc_label, collapse = ", "), ") were found in ",
        "metadata$", col_group, ". No LOESS trend line will be drawn."
      )
      col_qc_label <- NULL
    }
  }

  if (!is.null(arrange_group)) {
    if (!is.character(arrange_group)) {
      stop("'arrange_group' must be a character vector.")
    }
  }

  if (!is.null(plot_title)) {
    if (!is.character(plot_title) || length(plot_title) != 1L) {
      stop("'plot_title' must be a single character string or NULL.")
    }
  }

  if (identical(x, y)) {
    warning("'x' and 'y' are identical. No transformation may have been applied.")
  }

  # Seed validation (actual seeding done via .with_seed where needed)
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed))
      stop("'seed' must be a single numeric value or NULL.")
  }

  # -------------------------------------------------------------------------
  # Font size resolution
  # Priority: explicit arg > global_font_size scale > base_size-derived default
  # -------------------------------------------------------------------------
  gfs_factor      <- if (!is.null(global_font_size)) global_font_size / 11 else NULL
  r_title_size    <- .resolve_size(title_size,    gfs_factor, base_size + 1)
  r_subtitle_size <- .resolve_size(subtitle_size, gfs_factor, base_size + 2)
  r_xlab_size     <- .resolve_size(xlab_size,     gfs_factor, base_size)
  r_ylab_size     <- .resolve_size(ylab_size,     gfs_factor, base_size)
  r_axis_text     <- .resolve_size(NULL,          gfs_factor, base_size - 1)
  r_legend_title  <- .resolve_size(NULL,          gfs_factor, base_size)
  r_legend_text   <- .resolve_size(NULL,          gfs_factor, base_size - 1)
  r_strip_text    <- .resolve_size(NULL,          gfs_factor, base_size)

  # -------------------------------------------------------------------------
  # Base theme construction (mirrors plot_score approach)
  # -------------------------------------------------------------------------
  base_theme <- switch(
    theme_choice,
    "nature"  = ggplot2::theme_bw(base_size = base_size, base_family = font_family) +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_line(color = "grey92", linewidth = 0.4),
        panel.grid.minor = ggplot2::element_blank(),
        panel.border     = ggplot2::element_rect(color = "grey70", fill = NA, linewidth = 0.6),
        strip.background = ggplot2::element_rect(fill = "grey95", color = "grey70"),
        strip.text       = ggplot2::element_text(face = "bold", family = font_family,
                                                 size = r_strip_text)
      ),
    "minimal" = ggplot2::theme_minimal(base_size = base_size, base_family = font_family),
    "classic" = ggplot2::theme_classic(base_size = base_size, base_family = font_family),
    "bw"      = ggplot2::theme_bw(base_size     = base_size, base_family = font_family),
    "light"   = ggplot2::theme_light(base_size  = base_size, base_family = font_family),
    "dark"    = ggplot2::theme_dark(base_size   = base_size, base_family = font_family)
  )

  custom_theme <- ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    plot.title       = ggplot2::element_text(size = r_title_size,   hjust = 0.5,
                                             family = font_family),
    axis.title.x     = ggplot2::element_text(size = r_xlab_size,    family = font_family),
    axis.title.y     = ggplot2::element_text(size = r_ylab_size,    family = font_family),
    axis.text        = ggplot2::element_text(size = r_axis_text,    family = font_family),
    legend.title     = ggplot2::element_text(size = r_legend_title, family = font_family),
    legend.text      = ggplot2::element_text(size = r_legend_text,  family = font_family),
    strip.text       = ggplot2::element_text(size = r_strip_text,   face = "bold",
                                             family = font_family)
  )

  # -------------------------------------------------------------------------
  # Derive group palette dynamically (Okabe-Ito, colorblind-safe)
  # -------------------------------------------------------------------------
  if (is.null(arrange_group)) {
    all_groups_present <- sort(unique(as.character(metadata[[col_group]])))
    qc_p <- intersect(col_qc_label, all_groups_present)
    other_p <- setdiff(all_groups_present, qc_p)
    group_levels <- c(qc_p, other_p)
  } else {
    group_levels <- arrange_group
  }

  # Enforce factor levels in metadata for downstream plotting
  metadata[[col_group]] <- factor(metadata[[col_group]], levels = group_levels)

  okabe_ito <- c(
    "#0072B2", "#D55E00", "#009E73", "#56B4E9",
    "#E69F00", "#CC79A7", "#F0E442", "#000000"
  )
  n_groups     <- length(group_levels)
  group_colors <- stats::setNames(
    okabe_ito[seq_len(min(n_groups, length(okabe_ito)))],
    group_levels[seq_len(min(n_groups, length(okabe_ito)))]
  )
  if (n_groups > length(okabe_ito)) {
    extra_colors <- grDevices::hcl.colors(
      n_groups - length(okabe_ito), palette = "Dark 3"
    )
    group_colors <- c(
      group_colors,
      stats::setNames(extra_colors,
                      group_levels[(length(okabe_ito) + 1L):n_groups])
    )
  }

  # -------------------------------------------------------------------------
  # Derive batch shape palette (standard R shapes, up to 20 unique)
  # -------------------------------------------------------------------------
  batch_vec    <- as.factor(metadata[[col_batch]])
  batch_levels <- levels(batch_vec)

  # Prioritize solid shapes for visibility, then symbols
  shape_pool <- c(16, 17, 15, 18, 3, 4, 8, 1, 2, 0, 5, 6, 7, 10, 11, 12, 13, 14, 9, 20)
  n_batches_total <- length(batch_levels)

  if (n_batches_total > length(shape_pool)) {
    batch_shapes <- rep(shape_pool, length.out = n_batches_total)
  } else {
    batch_shapes <- shape_pool[seq_len(n_batches_total)]
  }
  batch_shapes <- stats::setNames(batch_shapes, batch_levels)

  # -------------------------------------------------------------------------
  # Build per-feature panels
  # -------------------------------------------------------------------------
  n_rows        <- nrow(x)
  injection_seq <- metadata[[col_injection]]
  # batch_vec moved above for palette derivation
  group_vec     <- metadata[[col_group]]
  y_axis_label  <- if (log_transform) "log10(Value)" else "Value"

  panel_labels        <- c(x_label, y_label)
  names(panel_labels) <- c("before", "after")

  plot_list <- lapply(plot_what, function(feat) {

    before_vals <- x[[feat]]
    after_vals  <- y[[feat]]

    if (log_transform) {
      before_vals <- log10(pmax(before_vals, 1e-10))
      after_vals  <- log10(pmax(after_vals,  1e-10))
    }

    plot_data <- data.frame(
      injection = rep(injection_seq, 2L),
      value     = c(before_vals, after_vals),
      panel     = rep(c("before", "after"), each = n_rows),
      group     = rep(group_vec,  2L),
      batch     = rep(batch_vec,  2L),
      stringsAsFactors = FALSE
    )
    plot_data$panel <- factor(plot_data$panel, levels = c("before", "after"))

    # LOESS eligibility check
    draw_loess <- !is.null(col_qc_label)
    if (draw_loess) {
      n_qc_obs <- sum(group_vec %in% col_qc_label)
      if (n_qc_obs < 3L) {
        warning(
          "Feature '", feat, "': fewer than 3 observations in the '",
          col_qc_label, "' group (", n_qc_obs, " found). ",
          "LOESS line suppressed for this feature."
        )
        draw_loess <- FALSE
      }
    }

    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = injection, y = value)
    ) +
      ggplot2::geom_point(
        ggplot2::aes(color = group, shape = batch),
        size  = point_size,
        alpha = point_alpha
      )

    if (draw_loess) {
      qc_data <- plot_data[plot_data$group %in% col_qc_label, , drop = FALSE]
      p <- p + ggplot2::geom_smooth(
        data      = qc_data,
        method    = NULL,
        formula   = NULL,
        se        = TRUE,
        na.rm     = TRUE,
        alpha     = 0.3,
        color     = "darkred",
        linetype  = "dashed",
        fullrange = TRUE
      )
    }

    p +
      ggplot2::facet_wrap(
        ~ panel,
        scales   = "free_y",
        labeller = ggplot2::labeller(panel = panel_labels)
      ) +
      ggplot2::scale_color_manual(
        values = group_colors,
        limits = group_levels,
        drop   = FALSE
      ) +
      ggplot2::scale_shape_manual(
        values = batch_shapes,
        limits = batch_levels,
        drop   = FALSE
      ) +
      ggplot2::labs(
        title = feat,
        x     = "Injection Sequence",
        y     = y_axis_label,
        color = "Group",
        shape = "Batch"
      ) +
      base_theme +
      custom_theme
  })

  # -------------------------------------------------------------------------
  # Main title and layout
  # -------------------------------------------------------------------------
  n_batches  <- length(unique(metadata[[col_batch]]))
  if (is.null(plot_title)) {
    main_title <- if (n_batches <= 1L) {
      paste0("Features: ", x_label, " vs ", y_label, " (Single Batch)")
    } else {
      paste0("Features: ", x_label, " vs ", y_label, " (", n_batches, " Batches)")
    }
  } else {
    main_title <- plot_title
  }

  n_plots <- length(plot_list)
  if (is.null(plot_cols)) {
    plot_cols <- if (n_plots == 1L) 1L else if (n_plots <= 6L) 2L else 3L
  }

  # Extract legend from the first plot and remove from all subplots to reduce clutter
  shared_legend <- NULL
  if (n_plots > 0) {
    tmp_grob <- ggplot2::ggplotGrob(plot_list[[1]])
    leg_idx  <- which(sapply(tmp_grob$grobs, function(x) x$name) == "guide-box")
    if (length(leg_idx) > 0) {
      shared_legend <- tmp_grob$grobs[[leg_idx]]
    }

    # Remove legends from all individual panels to maximize plot area
    plot_list <- lapply(plot_list, function(p) {
      p + ggplot2::theme(legend.position = "none")
    })
  }

  # Arrange subplots into a grid
  plots_grid <- gridExtra::arrangeGrob(grobs = plot_list, ncol = plot_cols)

  # Final arrangement with shared legend on the right
  main_title_grob <- grid::textGrob(
    main_title,
    gp = grid::gpar(
      fontsize   = r_subtitle_size,
      fontface   = "bold",
      fontfamily = font_family
    )
  )

  if (!is.null(shared_legend)) {
    combined_plot <- gridExtra::arrangeGrob(
      plots_grid, shared_legend,
      ncol   = 2,
      widths = c(10, 1.5), # Maximize plot space, narrow legend column
      top    = main_title_grob
    )
  } else {
    combined_plot <- gridExtra::arrangeGrob(
      plots_grid,
      top = main_title_grob
    )
  }

  grid::grid.newpage()
  grid::grid.draw(combined_plot)

  message("Before/after comparison plots created successfully.")
  return(invisible(combined_plot))
}