# Package-wide constants and R CMD check bindings.
# ponytail: one source of truth for method vocabularies, so run_DIpreprocess()
# pre-flight checks and the individual run_* validators can never drift apart.

#' @importFrom ggplot2 .data
NULL

# Supported method names (lower-case), shared by validators and docs.
.normalize_methods <- c("none", "sum", "median", "specific_factor",
                        "pqn_global", "pqn_reference", "pqn_group",
                        "pqn", "quantile",
                        "col_rel_abundance", "row_rel_abundance")
.transform_methods <- c("none", "log2", "log10", "sqrt", "cbrt", "clr",
                        "arcsin_sqrt", "vsn", "glog")
.scale_methods     <- c("none", "mean", "auto", "pareto", "ln",
                        "vast", "svast", "xvast")
.plot_themes       <- c("nature", "minimal", "classic", "bw", "light", "dark")
