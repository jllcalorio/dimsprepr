#' Plot Pairwise Before-and-After Comparison for Selected Features
#'
#' @description
#' Generates a multi-panel plot comparing the signal of selected features
#' before and after any data transformation, correction, or manipulation.
#' Each panel displays one feature across the injection sequence, colored by
#' sample group and shaped by batch, with an optional LOESS trend line fitted
#' to a user-defined reference group (e.g., QC samples).
#'
#' @param x A `data.frame`, `tibble`, or `matrix` representing the data
#'   \emph{before} the transformation. Rows are samples; columns are features.
#'   Column names may contain special characters. Alternatively, a
#'   [run_DIpreprocess()] result: `x` then becomes its missing-filtered data
#'   (`$data_missing_filtered`), `y` its drift/batch-corrected data
#'   (`$data_corrected`), and `metadata` its standardized metadata, so `y` and
#'   `metadata` can be omitted.
#' @param y A `data.frame`, `tibble`, or `matrix` representing the data
#'   \emph{after} the transformation. Must have the same dimensions and column
#'   names as `x`. Ignored when `x` is a `run_DIpreprocess` result.
#' @param metadata A `data.frame` containing sample-level metadata shared by
#'   both `x` and `y`. Must have the same number of rows as `x` and `y`.
#'   Optional when `x` is a `run_DIpreprocess` result.
#' @param plot_what Character vector or `NULL`. Column names of `x` to plot.
#'   If `NULL` (default), up to 4 features are selected at random. Features
#'   not found are skipped with a warning; if none is found, the function
#'   falls back to random selection.
#' @param col_injection Character. Column in `metadata` with the injection
#'   sequence (numeric run order). Default: `"InjectionSequence"` (always used
#'   for a `run_DIpreprocess` input).
#' @param col_batch Character. Column in `metadata` with batch identifiers.
#'   Default: `"Batch"` (always used for a `run_DIpreprocess` input).
#' @param col_group Character. Column in `metadata` with sample group labels.
#'   Default: `"Group"` (always used for a `run_DIpreprocess` input).
#' @param col_qc_label Character vector or `NULL`. Value(s) of `col_group`
#'   identifying the reference group for the LOESS trend line, typically QC
#'   samples. `NULL` draws no trend line. Default: `c("QC", "SQC", "EQC")`;
#'   for a `run_DIpreprocess` input, its `qc_types` are used.
#' @param arrange_group Character vector or `NULL`. Custom order for the
#'   sample groups in `col_group`. If `NULL` (default), groups are arranged
#'   alphanumerically, with `col_qc_label` values placed first.
#' @param plot_title Character or `NULL`. Main title for the combined plot.
#'   If `NULL` (default), a title is generated from `x_label`, `y_label`, and
#'   the number of batches.
#' @param x_label Character. Label of the before panel. Default: `"Before"`.
#' @param y_label Character. Label of the after panel. Default: `"After"`.
#' @param log_transform Logical. If `TRUE` (default), plot `log10` values
#'   (floored at `1e-10`).
#' @param point_size Numeric. Point size. Default: `1.5`.
#' @param point_alpha Numeric in (0, 1\]. Point transparency. Default: `0.7`.
#' @param plot_cols Integer or `NULL`. Number of columns of the panel grid;
#'   `NULL` (default) chooses it from the number of features.
#' @param theme Character. ggplot2 theme for each panel: `"nature"`
#'   (default; a clean style based on `theme_bw()`), `"minimal"`,
#'   `"classic"`, `"bw"`, `"light"`, or `"dark"`.
#' @param base_size Numeric. Base font size in points. All text sizes follow
#'   it unless set explicitly below. Default: `11`, ggplot2's default, which
#'   suits the default 7 x 5 inch figures of R Markdown and Quarto.
#' @param font_family Character. Font family for all text. Default: `"sans"`.
#' @param seed Numeric or `NULL`. Seed for the random feature selection; the
#'   global random number state is restored afterwards. Default: `123`.
#' @param title_size Numeric. Font size of the feature panel titles.
#'   Default: `base_size`.
#' @param subtitle_size Numeric. Font size of the main title above all
#'   panels. Default: `base_size * 1.2`.
#' @param xlab_size,ylab_size Numeric. Font sizes of the x- and y-axis titles.
#'   Default: `base_size`.
#'
#' @details
#' `x` and `y` must share identical row counts, column counts, and column
#' names. Row order is assumed to correspond between `x`, `y`, and
#' `metadata`. Missing values are left out of the panels.
#'
#' The LOESS trend line (Cleveland, 1979) is fitted only to the group given in
#' `col_qc_label`; with fewer than 3 observations of that group a warning is
#' issued and the line is suppressed. Groups are coloured with the
#' colour-blind-safe Okabe-Ito palette (Okabe & Ito, 2008).
#'
#' Like a ggplot, the plot is drawn when it is printed, which happens
#' automatically when the call is not assigned, so it appears directly in the
#' console and in R Markdown and Quarto documents. Inside loops or functions,
#' call `print()` on it.
#'
#' @return The combined plot: a `gtable` (from [gridExtra::arrangeGrob()])
#'   with class `"dimsprepr_plot"`, drawn when printed; [ggplot2::ggsave()]
#'   can save it.
#'
#' @references
#' Cleveland, W.S. (1979). Robust locally weighted regression and smoothing
#' scatterplots. *Journal of the American Statistical Association*, 74(368),
#' 829-836. \doi{10.1080/01621459.1979.10481038}
#'
#' Okabe, M., & Ito, K. (2008). *Color Universal Design (CUD): How to make
#' figures and presentations that are friendly to colorblind people.*
#' \url{https://jfly.uni-koeln.de/color/}
#'
#' @family plotting functions
#' @seealso [run_driftBatchCorrect()], whose effect this plot is designed to
#'   show.
#'
#' @examples
#' data(dims_data)
#' data(dims_metadata)
#' raw <- replace(dims_data, dims_data == 0, NA)  # zeros = not detected
#' before <- raw[, 1:6]
#' after  <- run_driftBatchCorrect(before, dims_metadata,
#'                                 spline_smooth_param = 0.5, verbose = FALSE)$data
#'
#' plot_beforeafter(before, after, dims_metadata, plot_what = "mz_0001",
#'                  x_label = "Raw", y_label = "QCRSC-corrected")
#'
#' # A run_DIpreprocess() result can be passed directly as `x`; it then shows
#' # the data before and after drift correction (see the package vignette).
#'
#' @author John Lennon L. Calorio
#'
#' @export
plot_beforeafter <- function(
    x,
    y,
    metadata,
    plot_what     = NULL,
    col_injection = "InjectionSequence",
    col_batch     = "Batch",
    col_group     = "Group",
    col_qc_label  = c("QC", "SQC", "EQC"),
    arrange_group = NULL,
    plot_title    = NULL,
    x_label       = "Before",
    y_label       = "After",
    log_transform = TRUE,
    point_size    = 1.5,
    point_alpha   = 0.7,
    plot_cols     = NULL,
    theme         = "nature",
    base_size     = 11,
    font_family   = "sans",
    seed          = 123,
    title_size    = base_size,
    subtitle_size = base_size * 1.2,
    xlab_size     = base_size,
    ylab_size     = base_size
) {

  # --- Integration with run_DIpreprocess: drift-correction diagnostic ---
  if (inherits(x, "run_DIpreprocess")) {
    if (!is.null(x$error))
      stop("The run_DIpreprocess result contains an error: ", x$error, call. = FALSE)
    if (missing(metadata) || is.null(metadata)) metadata <- x$metadata
    col_qc_label <- if (is.null(x$parameters$qc_types)) c("QC", "SQC", "EQC") else
      x$parameters$qc_types
    # Missing-filtered (before) vs drift/batch-corrected (after), shared features
    common <- intersect(colnames(x$data_missing_filtered), colnames(x$data_corrected))
    y <- x$data_corrected[, common, drop = FALSE]
    x <- x$data_missing_filtered[, common, drop = FALSE]
    # Standardized column names of run_DIpreprocess()$metadata
    col_injection <- "InjectionSequence"
    col_batch     <- "Batch"
    col_group     <- "Group"
  }

  # -------------------------------------------------------------------------
  # Input validation
  # -------------------------------------------------------------------------
  .is_tabular <- function(z) is.data.frame(z) || is.matrix(z)
  if (!.is_tabular(x)) stop("'x' must be a data.frame, tibble, or matrix.", call. = FALSE)
  if (!.is_tabular(y)) stop("'y' must be a data.frame, tibble, or matrix.", call. = FALSE)

  if (!identical(dim(x), dim(y)))
    stop("'x' and 'y' must have the same dimensions. Got: x = ", nrow(x), " x ",
         ncol(x), ", y = ", nrow(y), " x ", ncol(y), ".", call. = FALSE)

  x_cols <- colnames(x)
  if (is.null(x_cols)) stop("'x' and 'y' must have named columns.", call. = FALSE)
  if (!identical(x_cols, colnames(y)))
    stop("'x' and 'y' must have identical column names in the same order.", call. = FALSE)

  x <- as.data.frame(x, check.names = FALSE)
  y <- as.data.frame(y, check.names = FALSE)

  if (!is.data.frame(metadata)) stop("'metadata' must be a data.frame.", call. = FALSE)
  if (nrow(metadata) != nrow(x))
    stop("'metadata' must have the same number of rows as 'x' and 'y'. Got: metadata = ",
         nrow(metadata), " rows, x = ", nrow(x), " rows.", call. = FALSE)

  for (col_arg in list(
    list(val = col_injection, name = "col_injection"),
    list(val = col_batch,     name = "col_batch"),
    list(val = col_group,     name = "col_group")
  )) {
    if (!is.character(col_arg$val) || length(col_arg$val) != 1L)
      stop("'", col_arg$name, "' must be a single character string.", call. = FALSE)
    if (!col_arg$val %in% names(metadata))
      stop("Column '", col_arg$val, "' not found in 'metadata'. Available columns: ",
           paste(names(metadata), collapse = ", "), ".", call. = FALSE)
  }

  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L || is.na(seed)))
    stop("'seed' must be a single numeric value or NULL.", call. = FALSE)

  plot_what <- .with_seed(seed, .resolve_plot_what(plot_what, x_cols))

  for (nm in c("point_size", "point_alpha", "base_size", "title_size",
               "subtitle_size", "xlab_size", "ylab_size"))
    .check_positive_numeric(get(nm), nm)
  if (point_alpha > 1) {
    warning("'point_alpha' is > 1; it will be clamped to 1.", call. = FALSE)
    point_alpha <- 1
  }

  if (!is.null(plot_cols)) {
    if (!is.numeric(plot_cols) || length(plot_cols) != 1L ||
        is.na(plot_cols) || plot_cols < 1L || plot_cols != floor(plot_cols))
      stop("'plot_cols' must be a single positive integer or NULL.", call. = FALSE)
    plot_cols <- as.integer(plot_cols)
  }

  if (!is.logical(log_transform) || length(log_transform) != 1L || is.na(log_transform))
    stop("'log_transform' must be a single logical value (TRUE or FALSE).", call. = FALSE)
  if (!is.character(font_family) || length(font_family) != 1L)
    stop("'font_family' must be a single character string.", call. = FALSE)

  if (!is.null(col_qc_label)) {
    if (!is.character(col_qc_label) || length(col_qc_label) < 1L)
      stop("'col_qc_label' must be a character vector or NULL.", call. = FALSE)
    if (!any(col_qc_label %in% metadata[[col_group]])) {
      warning("None of the 'col_qc_label' values (", paste(col_qc_label, collapse = ", "),
              ") were found in metadata$", col_group, ". No LOESS trend line will be drawn.",
              call. = FALSE)
      col_qc_label <- NULL
    }
  }

  if (!is.null(arrange_group) && !is.character(arrange_group))
    stop("'arrange_group' must be a character vector.", call. = FALSE)
  if (!is.null(plot_title) && (!is.character(plot_title) || length(plot_title) != 1L))
    stop("'plot_title' must be a single character string or NULL.", call. = FALSE)

  if (identical(x, y))
    warning("'x' and 'y' are identical. No transformation may have been applied.",
            call. = FALSE)

  base_theme   <- .build_base_theme(theme, base_size, font_family)
  custom_theme <- .build_custom_theme(title_size, xlab_size, ylab_size)

  # -------------------------------------------------------------------------
  # Group colours (Okabe-Ito, colour-blind safe) and batch shapes
  # -------------------------------------------------------------------------
  if (is.null(arrange_group)) {
    all_groups_present <- sort(unique(as.character(metadata[[col_group]])))
    qc_p <- intersect(col_qc_label, all_groups_present)
    group_levels <- c(qc_p, setdiff(all_groups_present, qc_p))
  } else {
    group_levels <- arrange_group
  }
  metadata[[col_group]] <- factor(metadata[[col_group]], levels = group_levels)

  okabe_ito <- c("#0072B2", "#D55E00", "#009E73", "#56B4E9",
                 "#E69F00", "#CC79A7", "#F0E442", "#000000")
  n_groups <- length(group_levels)
  group_colors <- stats::setNames(
    c(okabe_ito, grDevices::hcl.colors(max(0L, n_groups - length(okabe_ito)),
                                       palette = "Dark 3"))[seq_len(n_groups)],
    group_levels
  )

  batch_vec    <- as.factor(metadata[[col_batch]])
  batch_levels <- levels(batch_vec)
  # Solid shapes first for visibility, recycled beyond 20 batches
  shape_pool   <- c(16, 17, 15, 18, 3, 4, 8, 1, 2, 0, 5, 6, 7, 10, 11, 12, 13, 14, 9, 20)
  batch_shapes <- stats::setNames(rep(shape_pool, length.out = length(batch_levels)),
                                  batch_levels)

  # -------------------------------------------------------------------------
  # Per-feature panels
  # -------------------------------------------------------------------------
  n_rows        <- nrow(x)
  injection_seq <- metadata[[col_injection]]
  group_vec     <- metadata[[col_group]]
  y_axis_label  <- if (log_transform) "log10(Value)" else "Value"
  panel_labels  <- c(before = x_label, after = y_label)

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
      panel     = factor(rep(c("before", "after"), each = n_rows),
                         levels = c("before", "after")),
      group     = rep(group_vec, 2L),
      batch     = rep(batch_vec, 2L)
    )

    draw_loess <- !is.null(col_qc_label)
    if (draw_loess) {
      n_qc_obs <- sum(group_vec %in% col_qc_label & !is.na(before_vals))
      if (n_qc_obs < 3L) {
        warning("Feature '", feat, "': fewer than 3 observations in the '",
                paste(col_qc_label, collapse = ", "), "' group (", n_qc_obs,
                " found). LOESS line suppressed for this feature.", call. = FALSE)
        draw_loess <- FALSE
      }
    }

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$injection, y = .data$value)) +
      ggplot2::geom_point(ggplot2::aes(color = .data$group, shape = .data$batch),
                          size = point_size, alpha = point_alpha, na.rm = TRUE)

    if (draw_loess) {
      p <- p + ggplot2::geom_smooth(
        data      = plot_data[plot_data$group %in% col_qc_label, , drop = FALSE],
        method    = "loess",
        formula   = y ~ x,
        se        = TRUE,
        na.rm     = TRUE,
        alpha     = 0.3,
        color     = "darkred",
        linetype  = "dashed",
        fullrange = TRUE
      )
    }

    p +
      ggplot2::facet_wrap(~ panel, scales = "free_y",
                          labeller = ggplot2::labeller(panel = panel_labels)) +
      ggplot2::scale_color_manual(values = group_colors, limits = group_levels, drop = FALSE) +
      ggplot2::scale_shape_manual(values = batch_shapes, limits = batch_levels, drop = FALSE) +
      ggplot2::labs(title = feat, x = "Injection Sequence", y = y_axis_label,
                    color = "Group", shape = "Batch") +
      base_theme +
      custom_theme
  })

  # -------------------------------------------------------------------------
  # Main title, shared legend, and layout
  # -------------------------------------------------------------------------
  n_batches <- length(unique(metadata[[col_batch]]))
  main_title <- if (!is.null(plot_title)) plot_title else
    paste0("Features: ", x_label, " vs ", y_label,
           if (n_batches <= 1L) " (Single Batch)" else paste0(" (", n_batches, " Batches)"))

  n_plots <- length(plot_list)
  if (is.null(plot_cols)) plot_cols <- if (n_plots <= 2L) 1L else 2L

  # One legend for all panels, taken from the first
  tmp_grob <- ggplot2::ggplotGrob(plot_list[[1]])
  leg_idx  <- which(vapply(tmp_grob$grobs, function(g) identical(g$name, "guide-box"),
                           logical(1L)))
  shared_legend <- if (length(leg_idx) > 0L) tmp_grob$grobs[[leg_idx[1]]] else NULL
  plot_list <- lapply(plot_list, function(p) p + ggplot2::theme(legend.position = "none"))

  plots_grid <- gridExtra::arrangeGrob(grobs = plot_list, ncol = plot_cols)
  title_grob <- grid::textGrob(
    main_title,
    gp = grid::gpar(fontsize = subtitle_size, fontface = "bold", fontfamily = font_family)
  )

  combined_plot <- if (!is.null(shared_legend)) {
    gridExtra::arrangeGrob(plots_grid, shared_legend, ncol = 2,
                           widths = grid::unit.c(grid::unit(1, "null"),
                                                 sum(shared_legend$widths)),
                           top = title_grob)
  } else {
    gridExtra::arrangeGrob(plots_grid, top = title_grob)
  }
  class(combined_plot) <- c("dimsprepr_plot", class(combined_plot))
  combined_plot
}
