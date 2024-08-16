################## Age Independent CD4 Back Calculation ##################
################### Plot results of the MCMC model ####################

here::i_am("r/plots_ai.R")

# libraries
library(here)
library(tidyverse)

load(here("data/postproc_ai.RData"))

# tables
postproc_ai$undiag_prev |>
    filter(year %% 1 == 0.75) |>
    mutate(year = year - 0.75) |>
    filter(year > 2004) |>
    pivot_wider(names_from = "year", values_from = "value")

postproc_ai$incidence_annual |>
    filter(year > 2004) |>
    pivot_wider(names_from = "year", values_from = "value")

# figures
undiag_prev_plt <- bind_rows(
    hiv_2022_rw$undiag_prev |> mutate(group = " No counterfactual"),
    hiv_2022_rw_cf1$undiag_prev |> mutate(group = "a. Incidence unaffected"),
    hiv_2022_rw_cf2$undiag_prev |> mutate(group = "b. Diagnosis probabilities unaffected")
) |>
    pivot_wider(names_from = "quartile") |>
    ggplot() +
    aes(year, y = `50%`, ymin = `2.5%`, ymax = `97.5%`, color = group, fill = group) +
    scale_y_continuous(limits = c(0, 9000), expand = expansion(mult = c(0, 0.05))) +
    scale_x_continuous(limits = c(2011, NA), breaks = c(2010, 2012, 2014, 2016, 2018, 2020)) +
    geom_line(linewidth = 0.8) +
    geom_ribbon(alpha = 0.2, linetype = "blank") +
    # data = hiv_2022_rw_cf1$undiag_prev |> mutate(group = "Counterfactual diagnoses") |>
    # pivot_wider(names_from = "quartile")) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom") +
    labs(x = "Year", y = "Undiagnosed HIV prevalence", title = "Estimated undiagnosed HIV prevalence", subtitle = "Lines and shaded areas indicate posterior median and 95% CrI", color = "", fill = "") +
    scale_color_manual(values = c("#F8766D", "#00BA38", "#619CFF")) +
    scale_fill_manual(values = c("#F8766D", "#00BA38", "#619CFF"))

ggsave(here("figures/undiag_prev_22.pdf"), width = 8, height = 6.5, device = cairo_pdf)
ggsave(here("figures/undiag_prev_22.png"), width = 8, height = 6.5)

incidence_plt <- bind_rows(
    hiv_2022_rw$incidence |> mutate(group = " No counterfactual"),
    hiv_2022_rw_cf1$incidence |> mutate(group = "a. Incidence unaffected"),
    hiv_2022_rw_cf2$incidence |> mutate(group = "b. Diagnosis probabilities unaffected")
) |>
    pivot_wider(names_from = "quartile") |>
    ggplot() +
    aes(quarter, y = `50%`, ymin = `2.5%`, ymax = `97.5%`, color = group, fill = group) +
    scale_y_continuous(limits = c(0, 1000), expand = expansion(mult = c(0, 0.05))) +
    scale_x_continuous(limits = c(2011, NA), breaks = c(2010, 2012, 2014, 2016, 2018, 2020)) +
    geom_line(linewidth = 0.8) +
    geom_ribbon(alpha = 0.2, linetype = "blank") +
    # data = hiv_2021_spl_fk$incidence |> mutate(group = "Counterfactual diagnoses") |>
    # pivot_wider(names_from = "quartile")) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom") +
    labs(x = "Year", y = "HIV incidence", title = "Estimated HIV incidence", subtitle = "Lines and shaded areas indicate posterior median and 95% CrI", color = "", fill = "") +
    scale_color_manual(values = c("#F8766D", "#00BA38", "#619CFF")) +
    scale_fill_manual(values = c("#F8766D", "#00BA38", "#619CFF"))

ggsave(here("figures/incidence_22.pdf"), width = 8, height = 6.5, device = cairo_pdf)
ggsave(here("figures/incidence_22.png"), width = 8, height = 6.5)

diag_prob_plt <- bind_rows(
    hiv_2022_rw$diag_prob |> mutate(group = " No counterfactual"),
    hiv_2022_rw_cf1$diag_prob |> mutate(group = "a. Incidence unaffected"),
    hiv_2022_rw_cf2$diag_prob |> mutate(group = "b. Diagnosis probabilities unaffected")
) |>
    pivot_wider(names_from = "quartile") |>
    mutate(strata = factor(strata, levels = c("CD4<200", "CD4 200-350", "CD4 350-500", "CD4>500"))) |>
    ggplot() +
    aes(quarter, y = `50%`, ymin = `2.5%`, ymax = `97.5%`, color = group, fill = group) +
    scale_y_continuous(labels = scales::percent, limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    scale_x_continuous(limits = c(2011, NA), breaks = c(2010, 2012, 2014, 2016, 2018, 2020)) +
    geom_line(linewidth = 0.8) +
    geom_ribbon(alpha = 0.2, linetype = "blank") +
    # data = hiv_2021_spl_fk$diag_prob |> mutate(group = "Counterfactual diagnoses") |> pivot_wider(names_from = "quartile") |>
    # mutate(strata = factor(strata, levels = c("CD4<200", "CD4 200-350", "CD4 350-500", "CD4>500")))) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom") +
    labs(x = "Year", y = "Diagnosis probability", title = "Estimated HIV diagnosis probabilities, by CD4 strata", subtitle = "Lines and shaded areas indicate posterior median and 95% CrI", color = "", fill = "") +
    scale_color_manual(values = c("#F8766D", "#00BA38", "#619CFF")) +
    scale_fill_manual(values = c("#F8766D", "#00BA38", "#619CFF")) +
    facet_wrap(~strata)

ggsave(here("figures/diag_prob_22.pdf"), width = 8, height = 6.5, device = cairo_pdf)
ggsave(here("figures/diag_prob_22.png"), width = 8, height = 6.5)

# generate a tibble containing the HIV diagnosis data
load(here("data/test_data_ai.RData"))
hiv_dx <- tibble(
    hiv = stan_data$HIV + stan_data$AIDS,
    year = seq(1995, 2022.75, by = 0.25)
)

diag_plt <- bind_rows(
    hiv_2022_rw$diagnoses |> mutate(group = " No counterfactual"),
    hiv_2022_rw_cf1$diagnoses |> mutate(group = "a. Incidence unaffected"),
    hiv_2022_rw_cf2$diagnoses |> mutate(group = "b. Diagnosis probabilities unaffected")
) |>
    pivot_wider(names_from = "quartile") |>
    ggplot() +
    aes(year, y = `50%`, color = group, fill = group) +
    scale_y_continuous(limits = c(0, 1000), expand = expansion(mult = c(0, 0.05))) +
    scale_x_continuous(limits = c(2011, NA), breaks = c(2010, 2012, 2014, 2016, 2018, 2020)) +
    geom_line(linewidth = 0.8) +
    geom_ribbon(alpha = 0.2, linetype = "blank", aes(ymin = `2.5%`, ymax = `97.5%`)) +
    # data = hiv_2021_rw$diagnoses |> mutate(group = "Counterfactual diagnoses") |> pivot_wider(names_from = "quartile")) +
    geom_point(data = hiv_dx, aes(x = year, y = hiv), color = "red", fill = "red") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom") +
    labs(x = "Year", y = "HIV diagnoses", title = "Estimated HIV diagnoses", subtitle = "Lines and shaded areas indicate posterior median and 95% CrI, points indicate observed\nHIV diagnoses", color = "", fill = "") +
    scale_color_manual(values = c("#F8766D", "#00BA38", "#619CFF")) +
    scale_fill_manual(values = c("#F8766D", "#00BA38", "#619CFF"))

ggsave(here("figures/hiv_diag_22.pdf"), width = 8, height = 6.5, device = cairo_pdf)
ggsave(here("figures/hiv_diag_22.png"), width = 8, height = 6.5)

undiag_prev_plt
incidence_plt
diag_prob_plt
diag_plt
