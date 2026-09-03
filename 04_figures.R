options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

source("00_config.R")

font_family <- "Arial"
ink <- "#233238"
grid_col <- "#DDE5E8"
teal <- "#007C83"
orange <- "#D97706"
blue <- "#3274A1"
purple <- "#7A5195"
red <- "#B64545"

theme_paper <- function(base_size = 9) {
  theme_minimal(base_size = base_size, base_family = font_family) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.15), colour = ink, margin = margin(b = 4)),
      plot.subtitle = element_text(size = rel(0.92), colour = "#526269", margin = margin(b = 8)),
      plot.caption = element_text(size = rel(0.78), colour = "#64747A", hjust = 0, margin = margin(t = 7)),
      axis.title = element_text(face = "bold", colour = ink),
      axis.text = element_text(colour = ink),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.3, colour = grid_col),
      legend.title = element_text(face = "bold"),
      legend.position = "bottom",
      plot.margin = margin(8, 10, 8, 8)
    )
}

export_figure <- function(plot, stem, width, height) {
  ggsave(file.path(figures_dir, paste0(stem, ".svg")), plot, width = width, height = height,
         units = "in", device = svglite::svglite, bg = "white")
  ggsave(file.path(figures_dir, paste0(stem, ".pdf")), plot, width = width, height = height,
         units = "in", device = cairo_pdf, bg = "white")
  ragg::agg_tiff(file.path(figures_dir, paste0(stem, ".tiff")), width = width,
                 height = height, units = "in", res = 600, compression = "lzw")
  print(plot)
  dev.off()
  ragg::agg_png(file.path(figures_dir, paste0(stem, ".png")), width = width,
                height = height, units = "in", res = 220)
  print(plot)
  dev.off()
}

# Figure 1: cohort construction and event structure
flow <- data.table(
  step = 1:5,
  label = c(
    "Malignant C22.1/C24.0 records\n2000-2023",
    "Histology 8160/8162\n2000-2023",
    "Primary hazard cohort\n2004-2023",
    "Cause-known cohort\nfor cancer-specific analyses",
    "5-year absolute-risk cohort\n2004-2018"
  ),
  n = c(45243, 31086, 28654, 28323, 17751),
  x = c(1, 1, 1, 1, 2.25),
  y = c(5, 4, 3, 2, 3)
)
flow[, display := paste0(label, "\nN = ", comma(n))]
flow_export <- flow[, .(step, label = gsub("\\n", " ", label), n)]
fwrite(flow_export, file.path(source_dir, "figure1_cohort_flow.csv"))

edges <- data.table(
  x = c(1, 1, 1, 1.18), y = c(4.72, 3.72, 2.72, 3),
  xend = c(1, 1, 1, 2.02), yend = c(4.28, 3.28, 2.28, 3)
)
p_flow <- ggplot() +
  geom_segment(data = edges, aes(x, y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.10, "inches")), linewidth = 0.55, colour = "#77878D") +
  geom_label(data = flow, aes(x, y, label = display), size = 3.0, family = font_family,
             linewidth = 0.35, label.padding = unit(0.25, "lines"),
             fill = c("#EEF6F6", "#EEF6F6", "#DCEFF0", "#DCEFF0", "#FFF3E1"), colour = ink) +
  coord_cartesian(xlim = c(0.45, 2.85), ylim = c(1.55, 5.45), clip = "off") +
  theme_void(base_family = font_family) +
  labs(title = "A  Cohort construction") +
  theme(plot.title = element_text(face = "bold", colour = ink, size = 10.5), plot.margin = margin(8, 8, 8, 8))

events <- fread(file.path(tables_dir, "age_site_event_cells.csv"))
event_long <- melt(
  events,
  id.vars = c("age_group", "site_group", "n"),
  measure.vars = c("cancer_deaths", "other_deaths", "unknown_cause_deaths", "alive"),
  variable.name = "status", value.name = "count"
)
event_long[, proportion := count / n]
event_long[, group := factor(paste(age_group, site_group, sep = " / "),
  levels = c("15-39 / iCCA", "15-39 / eCCA", "40-64 / iCCA", "40-64 / eCCA", "65+ / iCCA", "65+ / eCCA"))]
event_long[, status := factor(status,
  levels = c("cancer_deaths", "other_deaths", "unknown_cause_deaths", "alive"),
  labels = c("Cancer death", "Other-cause death", "Unknown cause", "Alive at cutoff"))]
fwrite(event_long, file.path(source_dir, "figure1_event_structure.csv"))

p_events <- ggplot(event_long, aes(group, proportion, fill = status)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.2) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.03))) +
  scale_fill_manual(values = c("Cancer death" = red, "Other-cause death" = orange,
                               "Unknown cause" = "#AAB5B9", "Alive at cutoff" = teal)) +
  labs(title = "B  Outcome structure by age and site", x = NULL, y = "Proportion of cohort", fill = NULL,
       caption = "iCCA, intrahepatic; eCCA, extrahepatic cholangiocarcinoma.") +
  theme_paper(8.5) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(axis.text.x = element_text(angle = 32, hjust = 1), legend.position = "bottom")

fig1 <- p_flow + p_events + plot_layout(widths = c(1.05, 1.25))
export_figure(fig1, "Figure1_cohort_and_events", 7.4, 4.6)

# Figure 2: continuous age-specific relative effects
curves <- fread(file.path(source_dir, "figure2_relative_effect_curves.csv"))
curves <- curves[(outcome == "Cancer-specific hazard" & model == "Model 2") |
                 (outcome == "Overall mortality" & model == "Model 2") |
                 (outcome == "Cancer death subdistribution")]
curves[, panel := factor(fcase(
  outcome == "Cancer-specific hazard", "A  Cancer-specific hazard\nModel 2 (conditional)",
  outcome == "Overall mortality", "B  Overall mortality\nModel 2 (conditional)",
  default = "C  Cancer-death subdistribution\nUnadjusted Fine-Gray"
), levels = c("A  Cancer-specific hazard\nModel 2 (conditional)",
              "B  Overall mortality\nModel 2 (conditional)",
              "C  Cancer-death subdistribution\nUnadjusted Fine-Gray"))]

interaction <- fread(file.path(tables_dir, "interaction_tests.csv"))
panel_p <- c(
  "A  Cancer-specific hazard\nModel 2 (conditional)" = interaction[outcome == "Cancer-specific hazard" & model == "Model 2", p],
  "B  Overall mortality\nModel 2 (conditional)" = interaction[outcome == "Overall mortality" & model == "Model 2", p],
  "C  Cancer-death subdistribution\nUnadjusted Fine-Gray" = interaction[outcome == "Cancer death subdistribution", p]
)
p_labels <- data.table(panel = factor(names(panel_p), levels = levels(curves$panel)),
                       label = paste0("Interaction P ", format.pval(panel_p, digits = 2, eps = 1e-4)),
                       age = 17, estimate = 1.50)

p_relative <- ggplot(curves, aes(age, estimate)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "#77878D", linewidth = 0.45) +
  geom_ribbon(aes(ymin = lower95, ymax = upper95), fill = teal, alpha = 0.16) +
  geom_line(colour = teal, linewidth = 0.95) +
  geom_text(data = p_labels, aes(age, estimate, label = label), inherit.aes = FALSE,
            hjust = 0, vjust = 1, family = font_family, size = 2.7, colour = ink) +
  facet_wrap(~ panel, nrow = 1) +
  scale_x_continuous(breaks = c(20, 30, 40, 50, 65, 75, 85)) +
  scale_y_log10(breaks = c(0.3, 0.5, 0.7, 1, 1.5), labels = label_number(accuracy = 0.1)) +
  coord_cartesian(ylim = c(0.25, 1.55)) +
  labs(title = "Age continuously modifies the mortality contrast between anatomical sites",
       subtitle = "Relative effect for eCCA versus iCCA; shaded bands are 95% confidence intervals",
       x = "Age at diagnosis, years", y = "Relative effect (log scale)",
       caption = "Model 1 adjusts for demographics and diagnosis era. Model 2 additionally conditions on stage, grade, and treatment.\nThe Fine-Gray curve is unadjusted and addresses the cancer-death subdistribution.") +
  theme_paper(8.5) +
  theme(strip.text = element_text(face = "bold", colour = ink, lineheight = 1.05), legend.position = "none")
export_figure(p_relative, "Figure2_age_relative_effects", 7.5, 4.3)

# Figure 3: 5-year cancer-death cumulative incidence and absolute differences
cif <- fread(file.path(source_dir, "figure3_cif_curves.csv"))
cif[, site_group := factor(site_group, levels = c("iCCA", "eCCA"))]
p_cif <- ggplot(cif, aes(age, cif5, colour = site_group)) +
  geom_line(linewidth = 1.0) +
  scale_colour_manual(values = c("iCCA" = purple, "eCCA" = orange)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), breaks = seq(0.5, 0.9, 0.1)) +
  coord_cartesian(ylim = c(0.50, 0.95)) +
  scale_x_continuous(breaks = c(20, 30, 40, 50, 65, 75, 85)) +
  labs(title = "A  Model-based 5-year cumulative incidence",
       subtitle = "Unadjusted Fine-Gray point estimates, diagnoses 2004-2018",
       x = "Age at diagnosis, years", y = "5-year cancer-death cumulative incidence", colour = NULL) +
  theme_paper(8.5)

rd <- fread(file.path(tables_dir, "five_year_risk_differences_age_groups.csv"))
rd[, age_group := factor(age_group, levels = c("65+", "40-64", "15-39"))]
rd[, label := paste0(sprintf("%+.1f", 100 * risk_difference_eCCA_minus_iCCA), " pp")]
fwrite(rd, file.path(source_dir, "figure3_nonparametric_risk_differences.csv"))
p_rd <- ggplot(rd, aes(risk_difference_eCCA_minus_iCCA, age_group)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#77878D", linewidth = 0.45) +
  geom_errorbar(aes(xmin = rd_lower95, xmax = rd_upper95), orientation = "y", width = 0.12, linewidth = 0.65, colour = teal) +
  geom_point(size = 2.5, colour = teal) +
  geom_text(aes(label = label), nudge_y = 0.24, hjust = 0.5, family = font_family, size = 2.7, colour = ink) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(-0.30, 0.04), breaks = seq(-0.30, 0, 0.10)) +
  labs(title = "B  Nonparametric 5-year\nabsolute risk difference",
       subtitle = "eCCA minus iCCA; Aalen-Johansen estimates\nwith 95% confidence intervals",
       x = "Risk difference", y = "Age group") +
  theme_paper(8.5) + theme(legend.position = "none")

fig3 <- p_cif + p_rd + plot_layout(widths = c(1.08, 0.92)) +
  plot_annotation(caption = "Negative values indicate lower 5-year cancer-death incidence for eCCA. pp, percentage points.",
                  theme = theme(plot.caption = element_text(family = font_family, size = 7, colour = "#64747A", hjust = 0)))
export_figure(fig3, "Figure3_five_year_absolute_risk", 7.5, 4.2)

# Figure 4: cohort-definition and follow-up-time sensitivity analyses
sens <- fread(file.path(tables_dir, "sensitivity_analyses.csv"))
sens[, analysis := factor(analysis, levels = rev(analysis))]
sens[, p_label := paste0("Pint=", format.pval(interaction_p, digits = 2, eps = 1e-4))]
p_sens <- ggplot(sens, aes(hr_eCCA_vs_iCCA_age65, analysis)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "#77878D", linewidth = 0.45) +
  geom_errorbar(aes(xmin = lower95, xmax = upper95), orientation = "y", width = 0.12, colour = blue, linewidth = 0.65) +
  geom_point(colour = blue, size = 2.4) +
  geom_text(aes(x = 1.10, label = p_label), hjust = 0, family = font_family, size = 2.5, colour = "#526269") +
  scale_x_continuous(limits = c(0.78, 1.25), breaks = c(0.8, 0.9, 1.0, 1.1, 1.2)) +
  labs(title = "A  Cohort-definition sensitivity analyses", subtitle = "Adjusted cancer-specific HR at age 65",
       x = "HR for eCCA versus iCCA", y = NULL) +
  theme_paper(8.2) + theme(legend.position = "none")

pw <- fread(file.path(tables_dir, "piecewise_followup_effects.csv"))
pw[, followup_period := factor(followup_period, levels = c("0-12 months", "12-36 months", ">36 months"))]
pw_test <- fread(file.path(tables_dir, "piecewise_interaction_test.csv"))
p_pw <- ggplot(pw, aes(age, estimate, colour = followup_period, group = followup_period)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "#77878D", linewidth = 0.45) +
  geom_errorbar(aes(ymin = lower95, ymax = upper95), width = 1.0, alpha = 0.65, linewidth = 0.45) +
  geom_line(linewidth = 0.65) + geom_point(size = 1.8) +
  scale_colour_manual(values = c("0-12 months" = teal, "12-36 months" = orange, ">36 months" = purple),
                      labels = c("0-12", "12-36", ">36 months")) +
  scale_y_log10(breaks = c(0.2, 0.4, 0.7, 1, 1.5, 2.5), labels = label_number(accuracy = 0.1)) +
  coord_cartesian(ylim = c(0.15, 2.5)) +
  scale_x_continuous(breaks = c(20, 30, 40, 50, 65, 75)) +
  labs(title = "B  Follow-up-time sensitivity analysis",
       subtitle = paste0("Three-way age x site x time interaction P = ", sprintf("%.3f", pw_test$p)),
       x = "Age at diagnosis, years", y = "Adjusted cancer-specific HR (log scale)", colour = "Follow-up") +
  theme_paper(8.2) +
  guides(colour = guide_legend(title.position = "top", nrow = 1))

fig4 <- p_sens + p_pw + plot_layout(widths = c(1, 1.15)) +
  plot_annotation(caption = "All estimates compare eCCA with iCCA. Pint denotes the overall age-by-site interaction P value within each cohort definition.",
                  theme = theme(plot.caption = element_text(family = font_family, size = 7, colour = "#64747A", hjust = 0)))
export_figure(fig4, "Figure4_sensitivity_analyses", 7.5, 4.3)

cat("Figures exported to:", figures_dir, "\n")
