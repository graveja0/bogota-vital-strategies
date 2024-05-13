v_tot_qaly <-
  map2(l_m_M, payoffs,~({
    v_u_str <- .y$qalys[,v_names_states] %>% as.vector()
    t(.x[,v_names_states] %*% v_u_str) %*% (v_dwe * v_wcc)
  })) %>%
  unlist()

v_tot_cost <-
  map2(l_m_M, payoffs,~({
    v_c_str <- .y$costs[,v_names_states] %>% as.vector()
    t(.x[,v_names_states] %*% v_c_str) %*% (v_dwc * v_wcc)
  })) %>%
  unlist()

df_cea <- calculate_icers(cost       = v_tot_cost,
                          effect     = v_tot_qaly,
                          strategies = v_names_str)
## CEA table in proper format ----
table_cea <- format_table_cea(df_cea) # Function included in "R/Functions.R"; depends on the `scales` package
table_cea

## CEA frontier -----
#* Function included in "R/Functions.R"; depends on the `ggplot2`  and `ggrepel` packages.
#* The latest version can be found in `dampack` package
p_cea_frontier <-
  plot(df_cea, label = "all", txtsize = 16) +
  expand_limits(x = max(table_cea$QALYs) + 0.1) +
  theme(legend.position = c(0.8, 0.2))

df_cea
