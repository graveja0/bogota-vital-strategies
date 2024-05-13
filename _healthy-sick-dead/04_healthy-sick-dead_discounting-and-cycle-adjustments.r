## 1.3. Get discounting and within-cycle correction weights.

# 1.3.1 Discounting weights for Costs

# v_dwc  <- 1 / ((1 + (params$d_c )) ^ (0:params$n_cycles))
# v_dwe  <- 1 / ((1 + (params$d_e )) ^ (0:params$n_cycles))

v_dwc  <- 1 / ((1 + (params$d_e * params$n_cycle_length)) ^ (0:params$n_cycles))
v_dwe  <- 1 / ((1 + (params$d_c * params$n_cycle_length)) ^ (0:params$n_cycles))


# 1.3.2.Within-cycle correction (WCC) using Simpson's 1/3 rule

v_wcc <- gen_wcc(n_cycles = params$n_cycles,
                 method = "Simpson1/3")
