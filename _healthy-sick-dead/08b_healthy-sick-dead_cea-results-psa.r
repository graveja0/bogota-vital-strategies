v_tot_qaly_psa <-
  map2(l_m_M_psa, payoffs_psa, ~({
    map2(.x, .y,~({
      v_u_str <- .y$qalys[,v_names_states] %>% as.vector()
      t(.x[,v_names_states] %*% v_u_str) %*% (v_dwe * v_wcc)
    })) %>%
      unlist()
  }))

v_tot_cost_psa <-
  map2(l_m_M_psa, payoffs_psa, ~({
    map2(.x, .y,~({
      v_c_str <- .y$costs[,v_names_states] %>% as.vector()
      t(.x[,v_names_states] %*% v_c_str) %*% (v_dwe * v_wcc)
    })) %>%
      unlist()
  }))

df_cea_psa <-
  map2_progress(v_tot_cost_psa, v_tot_qaly_psa,~({
    calculate_icers(cost       = .x,
                    effect     = .y,
                    strategies = v_names_str)}))

##################
df_c <- v_tot_cost_psa %>% bind_rows %>% data.frame()
df_e <- v_tot_qaly_psa %>% bind_rows() %>% data.frame()
df_psa_input <- params_psa %>% bind_cols() %>% data.frame()
l_psa <- dampack::make_psa_obj(cost          = df_c,
                               effectiveness = df_e,
                               parameters    = df_psa_input,
                               strategies    = v_names_str)
l_psa$strategies <- v_names_str
colnames(l_psa$effectiveness)<- v_names_str
colnames(l_psa$cost)<- v_names_str

#* Vector with willingness-to-pay (WTP) thresholds.
v_wtp <- seq(0, 200000, by = 5000)

### Cost-Effectiveness Scatter plot ----
#* Function included in "R/Functions.R"; depends on `tidyr` and `ellipse` packages.
#* The latest version can be found in `dampack` package
plot.psa(l_psa) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  xlab("Effectiveness (QALYs)") +
  guides(col = guide_legend(nrow = 2)) +
  theme(legend.position = "bottom")


### Incremental cost-effectiveness ratios (ICERs) with probabilistic output ----
#* Compute expected costs and effects for each strategy from the PSA
#* Function included in "R/Functions.R". The latest version can be found in `dampack` package
df_out_ce_psa <- summary.psa(l_psa)

#* Function included in "R/Functions.R"; depends on the `dplyr` package
#* The latest version can be found in `dampack` package
df_cea_psa <- calculate_icers(cost       = df_out_ce_psa$meanCost,
                              effect     = df_out_ce_psa$meanEffect,
                              strategies = df_out_ce_psa$Strategy)
df_cea_psa

### Plot cost-effectiveness frontier with probabilistic output ----
#* Function included in "R/Functions.R"; depends on the `ggplot2`  and `ggrepel` packages.
#* The latest version can be found in `dampack` package
plot.icers(df_cea_psa)

## Cost-effectiveness acceptability curves (CEACs) and frontier (CEAF) ---
#* Functions included in "R/Functions.R". The latest versions can be found in `dampack` package
ceac_obj <- ceac(wtp = v_wtp, psa = l_psa)
#* Regions of highest probability of cost-effectiveness for each strategy
summary.ceac(ceac_obj)
#* CEAC & CEAF plot
plot.ceac(ceac_obj) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  theme(legend.position = c(0.82, 0.5))

## Expected Loss Curves (ELCs) ----
#* Function included in "R/Functions.R".The latest version can be found in `dampack` package
elc_obj <- calc_exp_loss(wtp = v_wtp, psa = l_psa)
elc_obj
#* ELC plot
plot.exp_loss(elc_obj, log_y = FALSE,
              txtsize = 16, xlim = c(0, NA), n_x_ticks = 14,
              col = "full") +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  # geom_point(aes(shape = as.name("Strategy"))) +
  scale_y_continuous("Expected Loss (Thousand $)",
                     breaks = dampack::number_ticks(10),
                     labels = function(x) x/1000) +
  theme(legend.position = c(0.4, 0.7))

## Expected value of perfect information (EVPI) ----
#* Function included in "R/Functions.R". The latest version can be found in `dampack` package
evpi <- dampack::calc_evpi(wtp = v_wtp, psa = l_psa)
#* EVPI plot
plot.evpi(evpi, effect_units = "QALY")
