# Final submission figure mapping

The renderer retains stable analysis-era filenames. Use this map for the finalized article:

| Final label | Content | Renderer selection | Output stem in `results/figures/` |
|---|---|---|---|
| Figure 1 | Cohort selection | flow | `figure_1_inclusion_exclusion_flow` |
| Figure 2 | Static/class mortality forest | forest | `figure_5_landmark_mortality_forest` |
| Figure 3 | Joint trajectories | trajectories | `figure_2_joint_trajectories` |
| Figure 4 | KM and ICU-death CIF | survival | `figure_6_day7_landmark_survival` |
| Figure 5 | Multisystem phenotypes | phenotypes | `figure_4_class_clinical_profiles` |
| Figure 6 | Clinical flows | clinical-flows | `figure_3_class_treatment_outcomes` |
| Figure 7 | Day-1 classifier | prediction | `figure_7_day1_class_prediction` |
| Figure S1 | K selection | supporting | `supp_figure_s1_gbmtm_k_selection` |
| Figure S2 | Subgroups | subgroups | `supp_figure_s2_day7_subgroups` |
| Figure S3 | Correlations | phenotypes | `supp_figure_s3_system_associations` |
| Figure S4 | First-day management | management | `supp_figure_s6_day1_management_profiles` |
| Figure S5 | Six-hour management | management-profiles | `supp_figure_s4_six_hour_management_profiles` |
| Figure S6 | First-day treatment associations within fixed trajectory classes | treatment | `supp_figure_s5_treatment_associations` |
| Figure S7 | First-day treatment associations across predicted trajectory probabilities | treatment-probability | `supp_figure_s7_treatment_probability` |

Each stem has `.pdf`, `.svg`, `.tif`, and `.png` exports. `clinical-flows` also emits the disposition diagnostic; final Figure 6 is the stem listed above.

## Exploratory treatment figure legends

**Figure S6. First-day treatment associations within fixed trajectory classes.** Class-specific hazard ratios (HRs) for 28-day mortality and Fine–Gray subdistribution HRs (sHRs) for ICU mortality, with 95% confidence intervals (CIs), are shown from the 24-hour landmark in MIMIC-IV (A–B) and MIMIC-III (C–D). The four first-day contrasts are any recorded norepinephrine-equivalent dose (NEE), maximum NEE >0.30 versus 0 < NEE ≤0.30 µg/kg/min among records with a quantified positive dose, net fluid balance >35 versus ≤35 mL/kg (post hoc), and recorded continuous renal replacement therapy (CRRT). Models adjust for age, sex, race/ethnicity and the Oxford Acute Severity of Illness Score without mechanical-ventilation points. No MIMIC-III Class 1 patient had recorded CRRT; the interaction test compares Classes 2 and 3 for that contrast. Arrowheads indicate intervals extending beyond the plotting limits. P values are nominal. Source: [fixed-class results](../data/aggregate/display/figure_s6_fixed.csv).

**Figure S7. First-day treatment associations across predicted trajectory probabilities.** Conditional HRs for 28-day mortality and sHRs for ICU mortality from the 24-hour landmark. (A) Any recorded versus no recorded NEE. (B) Maximum NEE >0.30 versus 0 < NEE ≤0.30 µg/kg/min among records with a quantified positive dose. (C) Net fluid balance >35 versus ≤35 mL/kg (post hoc). (D) Recorded versus no recorded CRRT. Green solid lines denote Class 2 and red dashed lines denote Class 3; shading shows pointwise 95% CIs. Each curve varies one predicted probability while holding the other at zero and is shown over the overlap of the exposure arms’ marginal 2.5th–97.5th percentile ranges. The strips show observed binned probability distributions; within each class, the lower-exposure or no-record group appears above the higher-exposure or recorded group. The interaction P value tests both treatment-by-probability slopes jointly (2 df); the slope-difference P value tests their equality (1 df). Both are nominal. Sources: [curves](../data/aggregate/display/figure_s7_four_treatment_curves.csv), [tests](../data/aggregate/display/figure_s7_four_treatment_tests.csv), [support bins](../data/aggregate/display/figure_s7_four_treatment_histograms.csv).

## Final supplementary table sources

Paths below are relative to `data/aggregate/`. These are analysis/source CSVs rather than typeset manuscript tables; joined endpoints remain separate where applicable.

| Final table | Content | Source CSVs |
|---|---|---|
| S1 | Parent baseline | `display/publication_parent_baseline.csv` |
| S2 | Calendar observation | `endpoints/all_icu_calendar_day_observation_counts.csv` |
| S3 | Missingness | `display/baseline_availability.csv`; `endpoints/all_icu_supp_table_baseline_missingness.csv` |
| S4 | Static outcomes | `endpoints/all_icu_cox_static_effects.csv`; `icu/fine_gray_estimates.csv` |
| S5 | K candidates | `trajectory/all_icu_gbmtm_candidate_metrics_k1_k6.csv` |
| S6 | Degree/pruning | `trajectory/all_icu_gbmtm_parameter_sensitivity_summary.csv`; `trajectory/all_icu_gbmtm_parameter_sensitivity_all_candidates.csv` |
| S7 | Optimization | `trajectory/all_icu_gbmtm_k3_multistart_stability.csv` |
| S8 | Class outcomes | `endpoints/all_icu_cox_gbmtm_trajectory_effects.csv`; `icu/fine_gray_estimates.csv` |
| S9 | Subgroup interactions | `endpoints/all_icu_cox_gbmtm_subgroup_interactions.csv`; `icu/subgroup_interaction_tests.csv` |
| S10 | Class baseline | `display/publication_class_baseline.csv` |
| S11 | Flows | `display/supp_disposition_summary.csv`; `display/supp_treatment_flow_summary.csv`; corresponding `_nodes.csv`/`_links.csv` |
| S12 | Overall classifier | `classifier/prediction_performance.csv` |
| S13 | Class-specific classifier | `classifier/class_specific_performance.csv` |
