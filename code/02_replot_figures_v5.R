## -----------------------------------------------------------------
## Re-plot manuscript simulation figures (v5): no Holm brackets, no n-labels, plain-number labels
## Uses the saved run 2025_09_02_results_sim_syntheticdata.RData
## (no re-simulation needed). Point RDATA to that file.
## Outputs: Fig_L2L3_bias.png, Fig_L1L3_bias.png, Fig_GDS4_bias.png
##          (+ *_muhat.png = estimated-score versions)
## -----------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(ggstatsplot)
})

RDATA  <- "2025_09_02_results_sim_syntheticdata.RData"
OUTDIR <- "."
options(scipen = 999)   # no scientific notation anywhere (kept for the whole session)
load(RDATA)   # resA_*, resB_*, resB4_* (only these objects are used below)

## ---- reshape simulation arrays to long format --------------------
to_metric_df <- function(res, metric) {
  tbl <- as.data.frame(as.table(res[, metric, , , drop = FALSE]))
  names(tbl) <- c("iter","metric_name","method","miss","value")
  tbl$metric_name <- NULL
  names(tbl)[names(tbl) == "value"] <- metric
  tbl
}
array_to_df <- function(res, scenario) {
  df <- Reduce(function(x, y) merge(x, y, by = c("iter","method","miss")),
               list(to_metric_df(res, "MuHat"), to_metric_df(res, "Bias"),
                    to_metric_df(res, "MuTrue"), to_metric_df(res, "UsableRate")))
  df$scenario <- scenario
  df
}
arrays_to_long <- function(lst) {
  lst <- Filter(Negate(is.null), lst)
  bind_rows(lapply(names(lst), function(s) array_to_df(lst[[s]], s)))
}

df_A <- arrays_to_long(list(MCAR = resA_MCAR, MAR = resA_MAR))
df_B <- arrays_to_long(list(MbD = resB_MbD, MCAR = resB_MCAR, MAR = resB_MAR))
df_C <- arrays_to_long(list(MCAR = resB4_MCAR$mean, MAR = resB4_MAR$mean))  # GDS-4: MCAR | MAR

## ---- plotting ----------------------------------------------------
library(ggrepel)

## method names as in the manuscript (GDS-4 suffix dropped; panel/axis already say GDS-4)
pretty <- c(CompleteCase="Complete Case", Rule8="Rule-8", MeanSub="Sample Mean Sub",
            IpsativeX15="Ipsative", MI="Multiple Imputation",
            CC4="Complete Case", IpsativeX4="Ipsative", MI4="Multiple Imputation")

prep <- function(df, scn_labels = NULL) {
  df$scenario <- factor(df$scenario, levels = intersect(c("MCAR","MAR","MbD"), unique(df$scenario)))
  if (!is.null(scn_labels)) levels(df$scenario) <- scn_labels[levels(df$scenario)]
  m <- as.character(df$method)
  df$method <- factor(pretty[m], levels = unique(pretty[m][order(match(m, names(pretty)))]), ordered = TRUE)
  df
}

## number format: fixed decimals, ASCII minus (safe in any font), never "-0.00", never scientific
fmt <- function(v, d = 2) {
  r <- round(v, d) + 0                      # +0 turns -0 into 0
  formatC(r, format = "f", digits = d)
}

## everything added to each ggstatsplot panel (d = decimals shown in the mean label)
components <- function(red0 = FALSE, d = 2) {
  c(
    if (red0) list(ggplot2::geom_hline(yintercept = 0, colour = "red", linewidth = 0.5)),
    list(
      ## own mean marker + label (replaces ggstatsplot's, which prints 3.62e-04 etc.)
      ggplot2::stat_summary(fun = mean, geom = "point", colour = "darkred", size = 5),
      ggrepel::geom_label_repel(
        stat = "summary", fun = mean, parse = TRUE,
        aes(label = paste0("widehat(mu)[mean]~'='~'", fmt(after_stat(y), d), "'")),
        size = 3.4, hjust = 0, nudge_x = 0.14, direction = "y", seed = 1,   # label starts right of the dot
        min.segment.length = Inf, label.size = 0.25, label.padding = 0.18,
        box.padding = 0.15),
      ggplot2::scale_y_continuous(labels = scales::label_number(accuracy = 0.01)),
      ggplot2::scale_x_discrete(labels = function(l) {
        l <- gsub("\\s*\\(n\\s*=\\s*[0-9]+\\)", "", gsub("\n", " ", l))   # drop "(n = 100)"
        l <- gsub("\\s*\\(GDS-4\\)", "", l)
        scales::label_wrap(11)(l)                                          # wrap long names
      }),
      ggplot2::labs(x = NULL)
    )
  )
}

fig <- function(df, y, ylab, scn_labels = NULL, red0 = FALSE, d = 2) {
  dd <- prep(df, scn_labels); dd <- dd[!is.na(dd[[y]]), ]
  grouped_ggbetweenstats(
    data = dd, x = method, y = !!rlang::sym(y), grouping.var = scenario,
    plotgrid.args = list(ncol = length(unique(dd$scenario))),
    results.subtitle = FALSE,
    pairwise.display = "none",                 # no Holm brackets
    sample.size.label = FALSE,
    centrality.plotting = FALSE,               # ggstatsplot's own mean label off (we draw ours)
    ggtheme = ggplot2::theme_bw(base_size = 12),
    xlab = NULL, ylab = ylab,
    ggplot.component = components(red0, d))
}
save <- function(p, f, w) ggsave(file.path(OUTDIR, f), p, width = w, height = 5.75, dpi = 300)

## Figure: L2-L3, GDS-15
save(fig(df_A, "Bias",  "Bias", red0 = TRUE), "Fig_L2L3_bias.png", 10)
save(fig(df_A, "MuHat", "GDS-15 score"),                                  "Fig_L2L3_muhat.png", 10)

## Figure: L1-L3, GDS-15 (panels titled '... + MbD')
lab_B <- c(MCAR = "MCAR + MbD", MAR = "MAR + MbD", MbD = "MbD")
save(fig(df_B, "Bias",  "Bias", lab_B, red0 = TRUE), "Fig_L1L3_bias.png", 12)
save(fig(df_B, "MuHat", "GDS-15 score", lab_B),                                 "Fig_L1L3_muhat.png", 12)

## Figure: L1-L3, GDS-4 (MCAR | MAR)
save(fig(df_C, "Bias",  "Bias", red0 = TRUE, d = 3), "Fig_GDS4_bias.png", 10)
save(fig(df_C, "MuHat", "Estimated GDS-4 score"),                                   "Fig_GDS4_muhat.png", 10)
