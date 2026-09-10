options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

source("00_config.R")

font_family <- "Arial"
ink <- "#203038"
muted <- "#617178"
grid_col <- "#DDE5E8"
teal <- "#007C83"
orange <- "#D97706"
blue <- "#3274A1"
purple <- "#7A5195"
red <- "#B64545"

theme_paper <- function(base_size = 8.3) {
  theme_minimal(base_size = base_size, base_family = font_family) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.12), colour = ink, margin = margin(b = 4)),
      plot.subtitle = element_text(size = rel(0.88), colour = muted, margin = margin(b = 7)),
      plot.caption = element_text(size = rel(0.76), colour = muted, hjust = 0, margin = margin(t = 7)),
      axis.title = element_text(face = "bold", colour = ink),
      axis.text = element_text(colour = ink),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.28, colour = grid_col),
      legend.title = element_text(face = "bold"),
      legend.position = "bottom",
      strip.text = element_text(face = "bold", colour = ink),
      plot.margin = margin(7, 8, 7, 7)
    )
}

export_figure <- function(plot, stem, width, height) {
  ggsave(
    file.path(figures_dir, paste0(stem, ".svg")), plot,
    width = width, height = height, units = "in",
    device = svglite::svglite, bg = "white"
  )
  ggsave(
    file.path(figures_dir, paste0(stem, ".pdf")), plot,
    width = width, height = height, units = "in",
    device = cairo_pdf, bg = "white"
  )
  ragg::agg_tiff(
    file.path(figures_dir, paste0(stem, ".tiff")),
    width = width, height = height, units = "in", res = 600,
    compression = "lzw"
  )
  print(plot)
  dev.off()
  ragg::agg_png(
    file.path(figures_dir, paste0(stem, ".png")),
    width = width, height = height, units = "in", res = 240
  )
  print(plot)
  dev.off()
}

exact_tests <- fread(file.path(tables_dir, "prepub_age85_interaction_tests.csv"))
era_tests <- fread(file.path(tables_dir, "prepub_era_threeway_tests.csv"))
schema_tests <- fread(file.path(tables_dir, "prepub_schema3_interaction_tests.csv"))
time_tests <- fread(file.path(tables_dir, "prepub_time_varying_tests.csv"))
nonlinear_tests <- fread(file.path(tables_dir, "enhancement_nonlinearity_binary.csv"))
full_period_tests <- fread(file.path(tables_dir, "enhancement_fullperiod_tests.csv"))

exact_p <- exact_tests[analysis == "Exact ages 15-84" & model == "Model 1", p]
era_p <- era_tests[horizon_months == 24 & model == "Model 1", p]
schema_p <- schema_tests[model == "Model 1", p]
time_p <- time_tests[model == "Model 1", p]
nonlinear_p <- nonlinear_tests[test == "Nonlinear interaction component" & model == "Model 1", p]
full_period_p <- full_period_tests[test == "Overall age-by-site interaction" & model == "Model 1", p]
format_p <- function(p) if (p < 0.001) "P < 0.001" else paste0("P = ", sprintf("%.3f", p))

# A: adjusted absolute risk differences under flexible and constant site-effect models.
rd <- fread(file.path(tables_dir, "enhancement_cif_constant_vs_flexible.csv"))
rd[, rd_percent := 100 * risk_difference_eCCA_minus_iCCA]
rd[, lower_percent := 100 * rd_lower95]
rd[, upper_percent := 100 * rd_upper95]
rd[, model_type := factor(model_type, levels = c("Age-specific spline interaction", "Constant site effect"))]

p_rd <- ggplot(rd, aes(age, rd_percent, colour = model_type, fill = model_type, group = model_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = muted, linewidth = 0.4) +
  geom_ribbon(aes(ymin = lower_percent, ymax = upper_percent), alpha = 0.10, colour = NA) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 1.9) +
  scale_colour_manual(values = c("Age-specific spline interaction" = teal, "Constant site effect" = muted)) +
  scale_fill_manual(values = c("Age-specific spline interaction" = teal, "Constant site effect" = muted)) +
  scale_x_continuous(breaks = sort(unique(rd$age))) +
  scale_y_continuous(labels = label_number(suffix = " pp", accuracy = 1)) +
  labs(
    title = "A  Constant versus age-specific absolute-risk estimates",
    subtitle = "Standardized 5-year cancer-death risk difference; bootstrap 95% CI",
    x = "Age at diagnosis, years", y = "eCCA minus iCCA", colour = NULL, fill = NULL
  ) +
  theme_paper() + guides(fill = "none", colour = guide_legend(nrow = 2))

# B: compare fixed-age effects in the primary and expanded observation periods.
primary <- fread(file.path(tables_dir, "prepub_exact_age_effects.csv"))[
  analysis == "Exact ages 15-84" & model == "Model 1",
  .(age, estimate, lower95, upper95)
]
primary[, period := "2004-2023 primary cohort"]
expanded <- fread(file.path(tables_dir, "enhancement_fullperiod_effects.csv"))[, .(age, estimate, lower95, upper95)]
expanded[, period := "2000-2023 expanded cohort"]
period_effects <- rbindlist(list(primary, expanded))
period_effects[, period := factor(period, levels = c("2004-2023 primary cohort", "2000-2023 expanded cohort"))]
fwrite(period_effects, file.path(source_dir, "figure5_period_comparison.csv"))
p_era <- ggplot(period_effects, aes(age, estimate, colour = period, group = period)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = muted, linewidth = 0.4) +
  geom_line(linewidth = 0.68) +
  geom_point(size = 1.7) +
  geom_errorbar(aes(ymin = lower95, ymax = upper95), width = 1.2, linewidth = 0.35, alpha = 0.65) +
  scale_colour_manual(values = c(teal, orange)) +
  scale_x_continuous(breaks = c(20, 30, 40, 50, 65, 75)) +
  scale_y_log10(breaks = c(0.4, 0.6, 0.8, 1.0, 1.3), labels = label_number(accuracy = 0.1)) +
  coord_cartesian(ylim = c(0.35, 1.35)) +
  labs(
    title = "B  Expanded-period sensitivity analysis",
    subtitle = paste0("2000-2023 overall interaction ", format_p(full_period_p)),
    x = "Age at diagnosis, years", y = "Cancer-specific HR (log scale)", colour = NULL
  ) +
  theme_paper() + guides(colour = guide_legend(nrow = 2, byrow = TRUE))

# C: reliable 2010+ EOD schema, separating perihilar and distal eCCA.
schema <- fread(file.path(source_dir, "figure5_schema3_curves.csv"))
schema[, contrast := factor(
  contrast,
  levels = c("Perihilar eCCA", "Distal eCCA"),
  labels = c("Perihilar eCCA vs iCCA", "Distal eCCA vs iCCA")
)]
p_schema <- ggplot(schema, aes(age, estimate, colour = contrast, fill = contrast)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = muted, linewidth = 0.4) +
  geom_ribbon(aes(ymin = lower95, ymax = upper95), alpha = 0.11, colour = NA) +
  geom_line(linewidth = 0.88) +
  scale_colour_manual(values = c("Perihilar eCCA vs iCCA" = orange, "Distal eCCA vs iCCA" = purple)) +
  scale_fill_manual(values = c("Perihilar eCCA vs iCCA" = orange, "Distal eCCA vs iCCA" = purple)) +
  scale_x_continuous(breaks = c(20, 40, 60, 75, 84)) +
  scale_y_log10(breaks = c(0.3, 0.5, 0.7, 1.0, 1.4), labels = label_number(accuracy = 0.1)) +
  coord_cartesian(ylim = c(0.25, 1.45)) +
  labs(
    title = "C  Three-category anatomic schema",
    subtitle = paste0("EOD 2010+; overall age-by-schema interaction ", format_p(schema_p)),
    x = "Age at diagnosis, years", y = "Cancer-specific HR (log scale)", colour = NULL, fill = NULL
  ) +
  theme_paper() + guides(colour = guide_legend(nrow = 2), fill = "none")

# D: time-varying age-by-site effects over six follow-up intervals.
tv <- fread(file.path(source_dir, "figure5_time_varying.csv"))
tv <- tv[age %in% c(30, 50, 65, 75)]
tv[, followup_period := factor(
  followup_period,
  levels = c("0-6", "6-12", "12-24", "24-36", "36-60", ">60")
)]
tv[, age_label := factor(
  paste0("Age ", age),
  levels = paste0("Age ", c(30, 50, 65, 75))
)]
p_tv <- ggplot(tv, aes(followup_period, estimate, colour = age_label, group = age_label)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = muted, linewidth = 0.4) +
  geom_line(linewidth = 0.62) +
  geom_point(size = 1.7) +
  geom_errorbar(aes(ymin = lower95, ymax = upper95), width = 0.12, linewidth = 0.38, alpha = 0.7) +
  scale_colour_manual(values = c(teal, orange, purple, blue)) +
  scale_y_log10(breaks = c(0.2, 0.4, 0.7, 1.0, 1.5, 3), labels = label_number(accuracy = 0.1)) +
  coord_cartesian(ylim = c(0.18, 3.2)) +
  labs(
    title = "D  Effect stability over follow-up",
    subtitle = paste0("Six-interval age-by-site-by-time interaction ", format_p(time_p)),
    x = "Follow-up, months", y = "Cancer-specific HR (log scale)", colour = NULL
  ) +
  theme_paper() + theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  guides(colour = guide_legend(nrow = 1))

figure5 <- (p_rd + p_era) / (p_schema + p_tv) +
  plot_annotation(
    title = "Age and detailed anatomy jointly define the prognostic site contrast",
    subtitle = paste0(
      "Primary exact-age (15-84 years) interaction ", format_p(exact_p),
      "; nonlinear component ", format_p(nonlinear_p)
    ),
    caption = paste0(
      "HRs compare eCCA with iCCA unless the schema-specific contrast is stated; negative risk differences indicate lower adjusted 5-year cancer-death risk for eCCA.\n",
      "Model 1 adjusts for demographics and diagnosis era; younger and late-follow-up estimates have wider confidence intervals."
    ),
    theme = theme(
      plot.title = element_text(family = font_family, face = "bold", size = 12, colour = ink),
      plot.subtitle = element_text(family = font_family, size = 8.5, colour = muted),
      plot.caption = element_text(family = font_family, size = 7, colour = muted, hjust = 0)
    )
  )

export_figure(figure5, "Figure5_prepublication_validity_checks", 7.6, 7.1)
cat("Figure 5 exported to:", figures_dir, "\n")
