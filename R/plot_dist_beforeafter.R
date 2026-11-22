#' Plot Distribution Comparison Before and After Data Transformation
#'
#' @description
#' Generates density and box plots comparing data distributions before and
#' after any transformation, normalization, scaling, or correction step.
#' Box plots can be drawn per sample or per feature. Returns the individual
#' `ggplot2` objects and the combined plot.
#'
#' @param x A `data.frame`, `tibble`, or `matrix` representing the data
#'   \emph{before} transformation. Rows are samples; columns are features.
#'   Column names may contain special characters. Alternatively, a
#'   [run_DIpreprocess()] result: `x` then becomes its imputed data
#'   (`$data_imputed`) and `y` its final NONPLS data (`$data_nonpls`), so `y`
#'   can be omitted.
#' @param y A `data.frame`, `tibble`, or `matrix` representing the data
#'   \emph{after} transformation. Must have the same dimensions and column
#'   names as `x`. Ignored when `x` is a `run_DIpreprocess` result.
#' @param group_by Character. Perspective for box plots: `"sample"` (default)
#'   draws one box per sample, `"feature"` one box per feature.
#' @param plot_what Character vector or `NULL`. Subset of features to include
#'   in all plots. `NULL` (default) uses all features.
#' @param n_random Integer or `NULL`. Maximum number of samples (or features,
#'   for `group_by = "feature"`) shown in the box plots, drawn at random.
#'   `NULL` shows all. Default: `20`.
#' @param x_label Character. Label of the before state. Default: `"Before"`.
#' @param y_label Character. Label of the after state. Default: `"After"`.
#' @param plot_title Character or `NULL`. Main title of the combined plot.
#'   `NULL` (default) generates one.
#' @param x_fill Character. Fill colour of the before plots. Default:
#'   `"#56B4E9"` (Okabe-Ito sky blue).
#' @param y_fill Character. Fill colour of the after plots. Default:
#'   `"#E69F00"` (Okabe-Ito orange).
#' @param point_alpha Numeric in (0, 1\]. Transparency of the density and box
#'   fills. Default: `0.6`.
#' @param theme Character. ggplot2 theme: `"nature"` (default), `"minimal"`,
#'   `"classic"`, `"bw"`, `"light"`, or `"dark"`.
#' @param base_size Numeric. Base font size in points. All text sizes follow
#'   it unless set explicitly below. Default: `11`, ggplot2's default, which
#'   suits the default 7 x 5 inch figures of R Markdown and Quarto.
#' @param font_family Character. Font family for all text. Default: `"sans"`.
#' @param plot_cols Integer or `NULL`. Number of columns of the combined plot.
#'   `NULL` (default) uses 2 (densities on top, box plots below).
#' @param seed Numeric or `NULL`. Seed for the random selection of
#'   `n_random`; the global random number state is restored afterwards.
#'   Default: `123`.
#' @param title_size Numeric. Font size of the panel titles.
#'   Default: `base_size`.
#' @param subtitle_size Numeric. Font size of the main title above all
#'   panels. Default: `base_size * 1.2`.
#' @param xlab_size,ylab_size Numeric. Font sizes of the x- and y-axis titles.
#'   Default: `base_size`.
#'
#' @details
#' `x` and `y` must share identical row counts, column counts, and column
#' names, with corresponding row order. Missing values are left out.
#'
#' When `plot_what` is given, only those features are used. `n_random`
#' subsamples the box plots only; the density plots always use all selected
#' data. Box plot items are ordered by their median in `x`.
#'
#' Like a ggplot, the combined plot is drawn when the result is printed,
#' which happens automatically when the call is not assigned, so it appears
#' directly in the console and in R Markdown and Quarto documents. Inside
#' loops or functions, call `print()` on it.
#'
#' @return A list of class `"dimsprepr_dist_plot"` containing:
#' \describe{
#'   \item{`plot_density_before`}{Density plot of `x`.}
#'   \item{`plot_density_after`}{Density plot of `y`.}
#'   \item{`plot_box_before`}{Box plot of `x`.}
#'   \item{`plot_box_after`}{Box plot of `y`.}
#'   \item{`plot_combined`}{A `gtable` of all four panels.}
#' }
#' Printing the list draws the combined plot.
#'
#' @examples
#' data(dims_data)
#' data(dims_metadata)
#' raw <- replace(dims_data, dims_data == 0, NA)  # zeros = not detected
#' before <- run_mvimpute(raw, verbose = FALSE)$data
#' after  <- run_normalize(before, dims_metadata, method = "pqn_group",
#'                         verbose = FALSE)$data
#'
#' # Sample-wise comparison on the log scale
#' res <- plot_dist_beforeafter(log10(before), log10(after),
#'                              x_label = "Raw", y_label = "PQN", n_random = 15)
#'
#' # Feature-wise, subset of features; panels are returned for re-use
#' res <- plot_dist_beforeafter(log10(before), log10(after), group_by = "feature",
#'                              plot_what = colnames(before)[1:10],
#'                              theme = "classic")
#' res$plot_box_after
#'
#' @author John Lennon L. Calorio
#'
#' @family plotting functions
#' @seealso [run_normalize()], [run_transform()], and [run_scale()], whose
#'   effects this plot is designed to show.
#'
#' @export
plot_dist_beforeafter <- function(
    x,
    y,
    group_by      = "sample",
    plot_what     = NULL,
    n_random      = 20L,
    x_label       = "Before",
    y_label       = "After",
    plot_title    = NULL,
    x_fill        = "#56B4E9",
    y_fill        = "#E69F00",
    point_alpha   = 0.6,
    theme         = "nature",
    base_size     = 11,
    font_family   = "sans",
    plot_cols     = NULL,
    seed          = 123,
    title_size    = base_size,
    subtitle_size = base_size * 1.2,
    xlab_size     = base_size,
    ylab_size     = base_size
) {

  # --- Integration with run_DIpreprocess: imputed vs final (NONPLS) data ---
  if (inherits(x, "run_DIpreprocess")) {
    if (!is.null(x$error))
      stop("The run_DIpreprocess result contains an error: ", x$error, call. = FALSE)
    common <- intersect(colnames(x$data_imputed), colnames(x$data_nonpls))
    y <- x$data_nonpls[, common, drop = FALSE]
    x <- x$data_imputed[, common, drop = FALSE]
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

  group_by <- tolower(group_by)
  if (length(group_by) != 1L || !group_by %in% c("sample", "feature"))
    stop("'group_by' must be \"sample\" or \"feature\".", call. = FALSE)

  if (!is.null(plot_what)) {
    if (!is.character(plot_what) || length(plot_what) < 1L)
      stop("'plot_what' must be a non-empty character vector or NULL.", call. = FALSE)
    missing_feats <- setdiff(plot_what, x_cols)
    if (length(missing_feats) > 0L)
      stop("The following feature(s) in 'plot_what' were not found in 'x'/'y': ",
           paste(missing_feats, collapse = ", "), ".", call. = FALSE)
    x <- x[, plot_what, drop = FALSE]
    y <- y[, plot_what, drop = FALSE]
  }

  for (nm in c("point_alpha", "base_size", "title_size", "subtitle_size",
               "xlab_size", "ylab_size"))
    .check_positive_numeric(get(nm), nm)
  if (point_alpha > 1) {
    warning("'point_alpha' is > 1; it will be clamped to 1.", call. = FALSE)
    point_alpha <- 1
  }

  if (!is.null(n_random)) {
    if (!is.numeric(n_random) || length(n_random) != 1L ||
        is.na(n_random) || n_random < 1L || n_random != floor(n_random))
      stop("'n_random' must be a single positive integer or NULL.", call. = FALSE)
    n_random <- as.integer(n_random)
  }

  if (!is.null(plot_cols)) {
    if (!is.numeric(plot_cols) || length(plot_cols) != 1L ||
        is.na(plot_cols) || plot_cols < 1L || plot_cols != floor(plot_cols))
      stop("'plot_cols' must be a single positive integer or NULL.", call. = FALSE)
    plot_cols <- as.integer(plot_cols)
  }

  if (!is.character(font_family) || length(font_family) != 1L)
    stop("'font_family' must be a single character string.", call. = FALSE)
  if (!is.null(plot_title) && (!is.character(plot_title) || length(plot_title) != 1L))
    stop("'plot_title' must be a single character string or NULL.", call. = FALSE)
  if (!is.character(x_fill) || length(x_fill) != 1L)
    stop("'x_fill' must be a single character color string.", call. = FALSE)
  if (!is.character(y_fill) || length(y_fill) != 1L)
    stop("'y_fill' must be a single character color string.", call. = FALSE)

  if (identical(x, y))
    warning("'x' and 'y' are identical. No transformation may have been applied.",
            call. = FALSE)

  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L || is.na(seed)))
    stop("'seed' must be a single numeric value or NULL.", call. = FALSE)

  base_theme   <- .build_base_theme(theme, base_size, font_family)
  custom_theme <- .build_custom_theme(title_size, xlab_size, ylab_size, title_face = "bold")

  # -------------------------------------------------------------------------
  # Long format (base R, no tidyr dependency)
  # -------------------------------------------------------------------------
  .to_long <- function(df) {
    n_rows <- nrow(df)
    ids    <- if (!is.null(rownames(df))) rownames(df) else as.character(seq_len(n_rows))
    data.frame(
      id    = rep(ids, times = ncol(df)),
      key   = rep(names(df), each = n_rows),
      value = unlist(df, use.names = FALSE),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }
  long_x <- .to_long(x)
  long_y <- .to_long(y)

  # -------------------------------------------------------------------------
  # Density plots (all selected features, all samples)
  # -------------------------------------------------------------------------
  .density <- function(long, fill, label)
    ggplot2::ggplot(long, ggplot2::aes(x = .data$value)) +
      ggplot2::geom_density(fill = fill, color = fill, alpha = point_alpha, na.rm = TRUE) +
      ggplot2::labs(title = paste("Distribution -", label), x = "Value", y = "Density") +
      base_theme + custom_theme

  p_density_before <- .density(long_x, x_fill, x_label)
  p_density_after  <- .density(long_y, y_fill, y_label)

  # -------------------------------------------------------------------------
  # Box plots (optionally subsampled)
  # -------------------------------------------------------------------------
  item_col    <- if (group_by == "sample") "id" else "key"
  box_x_label <- if (group_by == "sample") "Sample" else "Feature"
  all_items   <- unique(long_x[[item_col]])

  selected_items <- if (!is.null(n_random) && n_random < length(all_items)) {
    .with_seed(seed, sample(all_items, n_random))
  } else {
    all_items
  }
  n_shown <- length(selected_items)

  sub_x <- long_x[long_x[[item_col]] %in% selected_items, , drop = FALSE]
  sub_y <- long_y[long_y[[item_col]] %in% selected_items, , drop = FALSE]

  # Order items by median of 'before' data for consistency across panels
  med_order <- names(sort(tapply(sub_x$value, sub_x[[item_col]], stats::median, na.rm = TRUE)))
  med_order <- c(med_order, setdiff(selected_items, med_order))   # all-NA items last
  sub_x[[item_col]] <- factor(sub_x[[item_col]], levels = med_order)
  sub_y[[item_col]] <- factor(sub_y[[item_col]], levels = med_order)

  .box <- function(sub, fill, label)
    ggplot2::ggplot(sub, ggplot2::aes(x = .data[[item_col]], y = .data$value)) +
      ggplot2::geom_boxplot(fill = fill, color = "grey30", alpha = point_alpha,
                            outlier.size = 0.8, na.rm = TRUE) +
      ggplot2::coord_flip() +
      ggplot2::labs(title = paste0(box_x_label, " Profiles - ", label, " (n = ", n_shown, ")"),
                    x = box_x_label, y = "Value") +
      base_theme + custom_theme +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = ggplot2::rel(0.7)))

  p_box_before <- .box(sub_x, x_fill, x_label)
  p_box_after  <- .box(sub_y, y_fill, y_label)

  # -------------------------------------------------------------------------
  # Combined plot (redundant axis labels hidden)
  # -------------------------------------------------------------------------
  main_title <- if (!is.null(plot_title)) plot_title else
    paste0("Distribution Comparison: ", x_label, " vs ", y_label, " (",
           toupper(substring(group_by, 1, 1)), substring(group_by, 2), " Perspective)")

  no_x <- ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                         axis.text.x  = ggplot2::element_blank())
  no_y <- ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                         axis.text.y  = ggplot2::element_blank())

  combined_plot <- gridExtra::arrangeGrob(
    p_density_before + no_x, p_density_after + no_x + no_y,
    p_box_before,            p_box_after + no_y,
    ncol = if (!is.null(plot_cols)) plot_cols else 2L,
    top  = grid::textGrob(
      main_title,
      gp = grid::gpar(fontsize = subtitle_size, fontface = "bold", fontfamily = font_family)
    )
  )

  out <- structure(
    list(
      plot_density_before = p_density_before,
      plot_density_after  = p_density_after,
      plot_box_before     = p_box_before,
      plot_box_after      = p_box_after,
      plot_combined       = combined_plot
    ),
    class = c("dimsprepr_dist_plot", "dimsprepr_plot", "list")
  )
  out
}
