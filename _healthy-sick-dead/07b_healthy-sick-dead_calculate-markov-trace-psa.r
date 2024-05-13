# PSA
for (k in 1:n_sim) {
  for (t in 1:n_cycles) {
    res <-
      map2(l_m_M_psa[[k]],m_P_psa[[k]],~({
        .x[paste0(t-1),] %*% .y
      }))
    l_m_M_psa[[k]] <-
      map2(l_m_M_psa[[k]],res,~({
        .x[paste0(t),] <- .y
        .x
      }))
  }
}
