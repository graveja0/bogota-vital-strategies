
# Construct shell markov trace
l_m_M <-
  m_P %>%
  map(~({
    tmp <- .x[1,] %>% as.matrix() %>% t()
    tmp <- matrix(0,ncol=(ncol(tmp)),nrow=(n_cycles+1),dimnames = list(paste0(0:n_cycles),colnames(.x)))
    tmp[1,1] <- 1
    tmp
  }))

l_m_M_psa <- rep(list(l_m_M),n_sim)


# Fill out the markov trace
for (t in 1:n_cycles) {
  res <-
    map2(l_m_M,m_P,~({
      .x[paste0(t-1),] %*% .y
    }))
  l_m_M <-
    map2(l_m_M,res,~({
      .x[paste0(t),] <- .y
      .x
    }))
}


