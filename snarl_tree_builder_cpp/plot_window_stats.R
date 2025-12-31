#!/usr/bin/env Rscript

# Usage: Rscript plot_window_stats.R window_stats.tsv [-d lower_percentile] [-D upper_percentile] [output.pdf]
# 
# Generates one plot per top-level chain showing hap density 
# (#hap-informative leaf snarls / window length in bp) across windows.
#
# The plot helps identify appropriate thresholds for tagging 
# "hap-snarl dense" windows for chunk point selection.
#
# Arguments:
#   -d: Lower percentile for hap-dense classification (default: 50 = median)
#   -D: Upper percentile for hap-dense classification (default: 75 = 75th percentile)
#   Windows are labeled as hap-dense if their hap_density falls within the
#   percentile range calculated per chromosome.

library(ggplot2)
library(data.table)
library(scales)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript plot_window_stats.R window_stats.tsv [-d lower_percentile] [-D upper_percentile] [output.pdf]")
}

# Default values
lower_percentile <- 50.0
upper_percentile <- 75.0
tsv_path <- NULL
out_pdf <- "window_stats_plots.pdf"

# Parse arguments
i <- 1
while (i <= length(args)) {
  if (args[i] == "-d" && i + 1 <= length(args)) {
    lower_percentile <- as.numeric(args[i + 1])
    i <- i + 2
  } else if (args[i] == "-D" && i + 1 <= length(args)) {
    upper_percentile <- as.numeric(args[i + 1])
    i <- i + 2
  } else if (is.null(tsv_path)) {
    tsv_path <- args[i]
    i <- i + 1
  } else {
    out_pdf <- args[i]
    i <- i + 1
  }
}

if (is.null(tsv_path)) {
  stop("Error: window_stats.tsv file path is required")
}

cat("Lower percentile:", lower_percentile, "\n")
cat("Upper percentile:", upper_percentile, "\n")
cat("Input TSV:", tsv_path, "\n")
cat("Output PDF:", out_pdf, "\n")

# Read data
cat("Reading:", tsv_path, "\n")
dt <- fread(tsv_path)

# Check for expected columns - support both old and new formats
needed_old <- c("chain_id", "window_idx", "start_node", "end_node", 
                "leaf_snarl_count", "hap_informative_count", "window_length_bp", "path_name")
needed_new <- c(needed_old, "hap_density", "is_hap_dense")

# Check which format we have
if (all(needed_new %in% colnames(dt))) {
  cat("Detected new format with pre-computed hap_density and is_hap_dense\n")
  # Already has hap_density, just ensure it's numeric
  dt[, hap_density := as.numeric(hap_density)]
} else if (all(needed_old %in% colnames(dt))) {
  cat("Detected old format, computing hap_density\n")
  # Calculate hap density (#hap-informative leaf snarls / window length in bp)
  # Multiply by 1000 to get per-kb density for readability
  dt[, hap_density := (hap_informative_count / window_length_bp) * 1000]
  dt[, is_hap_dense := "no"]  # placeholder
} else {
  missing <- setdiff(needed_old, colnames(dt))
  stop(paste("Missing columns in TSV:", paste(missing, collapse=", ")))
}

# Handle cases where hap_density is 0 or NA
dt[is.na(hap_density) | !is.finite(hap_density), hap_density := 0]

# Calculate percentile-based thresholds per chromosome and label hap-dense windows
cat("\nCalculating percentile-based thresholds per chromosome...\n")
dt[, threshold_lower := quantile(hap_density[hap_density > 0], probs = lower_percentile / 100, na.rm = TRUE), by = path_name]
dt[, threshold_upper := quantile(hap_density[hap_density > 0], probs = upper_percentile / 100, na.rm = TRUE), by = path_name]

# Label windows as hap-dense if they fall within the percentile range
dt[, is_hap_dense := "no"]
dt[hap_density >= threshold_lower & hap_density <= threshold_upper, is_hap_dense := "yes"]

# For chromosomes with no non-zero hap_density, set thresholds to 0
dt[is.na(threshold_lower), threshold_lower := 0]
dt[is.na(threshold_upper), threshold_upper := 0]

cat("Labeled windows based on percentile range [", lower_percentile, ", ", upper_percentile, "]\n", sep="")
cat("Total hap-dense windows:", sum(dt$is_hap_dense == "yes"), "out of", nrow(dt), "\n")

# Calculate window positions along the chain
# Use start_node as position reference, or cumulative window lengths
dt[, window_length_bp := as.numeric(window_length_bp)]
dt <- dt[order(chain_id, window_idx)]
dt[, window_start_bp := cumsum(shift(window_length_bp, fill = 0)), by = chain_id]
dt[, window_mid_bp := window_start_bp + window_length_bp / 2]

# Filter to chains with at least some hap-informative snarls
chain_summary <- dt[, .(
  total_windows = .N,
  total_hap_informative = sum(hap_informative_count),
  total_leaf_snarls = sum(leaf_snarl_count),
  total_length_bp = sum(window_length_bp),
  mean_hap_density = mean(hap_density, na.rm = TRUE),
  max_hap_density = max(hap_density, na.rm = TRUE)
), by = .(chain_id, path_name)]

cat("\n=== Chain Summary ===\n")
print(chain_summary[order(-total_hap_informative)][1:min(20, .N)])

# Only plot chains that have some hap-informative snarls
chains_to_plot <- chain_summary[total_hap_informative > 0]$chain_id
dt_filtered <- dt[chain_id %in% chains_to_plot]

cat("\nPlotting", length(chains_to_plot), "chains with hap-informative snarls...\n")

if (nrow(dt_filtered) == 0) {
  cat("No chains with hap-informative snarls found. Exiting.\n")
  quit(status = 0)
}

# Calculate global density statistics for reference
global_density_stats <- dt_filtered[hap_density > 0, .(
  median = quantile(hap_density, 0.50),
  mean = mean(hap_density),
  lower_pct = quantile(hap_density, lower_percentile / 100),
  upper_pct = quantile(hap_density, upper_percentile / 100),
  q90 = quantile(hap_density, 0.90),
  q95 = quantile(hap_density, 0.95)
)]
cat("\n=== Hap Density Statistics (per kb, non-zero windows only) ===\n")
print(global_density_stats)
cat("\nUsing percentile range [", lower_percentile, ", ", upper_percentile, "] for hap-dense classification\n", sep="")

# Create PDF with two plots per chain (normal and log scale)
pdf(out_pdf, width = 12, height = 6)

for (cid in chains_to_plot) {
  sub <- dt_filtered[chain_id == cid][order(window_idx)]
  if (nrow(sub) == 0) next
  
  path_lab <- unique(sub$path_name)[1]
  total_hap <- sum(sub$hap_informative_count)
  total_leaf <- sum(sub$leaf_snarl_count)
  
  # Get percentile thresholds for this chromosome
  chain_threshold_lower <- unique(sub$threshold_lower)[1]
  chain_threshold_upper <- unique(sub$threshold_upper)[1]
  
  # Use is_hap_dense column for labeling
  sub[, is_dense := is_hap_dense == "yes"]
  color_labels <- c("Not hap-dense", "Hap-dense")
  
  # Helper function to create base plot
  create_base_plot <- function(sub, use_log = FALSE) {
    p <- ggplot(sub, aes(x = window_idx, y = hap_density)) +
      geom_line(color = "steelblue", linewidth = 0.8) +
      geom_point(aes(color = is_dense), size = 1, alpha = 0.7) +
      scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "firebrick"),
                         labels = color_labels,
                         name = "Hap Dense") +
      geom_hline(yintercept = chain_threshold_lower, linetype = "dashed", color = "orange", linewidth = 0.6) +
      geom_hline(yintercept = chain_threshold_upper, linetype = "dotted", color = "darkgreen", linewidth = 0.6)
    
    # Apply log transformation if requested
    if (use_log) {
      p <- p + scale_y_log10(labels = scales::comma_format())
    }
    
    scale_label <- if (use_log) " (log10 scale)" else ""
    
    p <- p +
      labs(
        title = paste0("Chain: ", cid, scale_label),
        subtitle = paste0("Path: ", path_lab, " | Windows: ", nrow(sub), 
                         " | Total hap-informative: ", total_hap, "/", total_leaf,
                         " | Hap-dense: ", sum(sub$is_dense)),
        x = "Window Index",
        y = paste0("Hap Density (hap-informative snarls per kb)", scale_label),
        caption = paste0("Orange dashed = ", lower_percentile, "th percentile (", round(chain_threshold_lower, 3), 
                        "), Green dotted = ", upper_percentile, "th percentile (", round(chain_threshold_upper, 3), ")")
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "gray40"),
        legend.position = "top",
        panel.grid.minor = element_blank()
      )
    
    return(p)
  }
  
  # Plot 1: Normal scale
  p1 <- create_base_plot(sub, use_log = FALSE)
  print(p1)
  
  # Plot 2: Log scale (only if there are non-zero values)
  if (any(sub$hap_density > 0)) {
    p2 <- create_base_plot(sub, use_log = TRUE)
    print(p2)
  }
}

dev.off()
cat("\nWrote plots to:", out_pdf, "\n")


cat("\n=== Hap-Dense Classification Summary ===\n")
cat("Percentile range used: [", lower_percentile, ", ", upper_percentile, "]\n", sep="")
cat("Global thresholds (across all chromosomes):\n")
cat("  - ", lower_percentile, "th percentile: ", round(global_density_stats$lower_pct, 4), " per kb\n", sep="")
cat("  - ", upper_percentile, "th percentile: ", round(global_density_stats$upper_pct, 4), " per kb\n", sep="")
cat("\nWindows are labeled as hap-dense if their density falls within the\n")
cat("percentile range calculated per chromosome.\n")
cat("\nTotal hap-dense windows: ", sum(dt$is_hap_dense == "yes"), " (", 
    round(100 * sum(dt$is_hap_dense == "yes") / nrow(dt), 1), "% of all windows)\n", sep="")
