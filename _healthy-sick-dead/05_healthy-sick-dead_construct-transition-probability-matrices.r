##############################################
### Step 2: Process Transition Probabilities
##############################################

# 2.0. Function to obtain transition probability matrices from
# underlying rate matrix structure (as specified in Excel)

get_transition_matrices <- function(strategies, params, jumpover = NULL) {
  m_P <-
    strategies %>% map(~({
      n <- params %>% map_dbl(~(length(.x))) %>% max()
      m_Q <- P_ <- list()

      m_Q_ <-
        transition_rates_raw %>%
        filter(strategy == .x) %>%
        filter(from %in% all_of(v_names_states))   %>%
        select_at(vars(from,v_names_states)) %>%
        column_to_rownames(var = "from") %>%
        as.matrix(); m_Q_

      if (!is.null(jumpover)) {
        for (jj in names(jumpover)) {
          ff <- jumpover[[jj]][[1]]; ff
          tt <- jumpover[[jj]][[2]]; tt
          tmp_ <- rbind(cbind(m_Q_,rep("0",nrow(m_Q_))),rep("0",ncol(m_Q_)+1)); tmp_
          colnames(tmp_) <- c(colnames(m_Q_),jj)
          rownames(tmp_) <- c(rownames(m_Q_),jj)
          tmp_[ff,jj] <- m_Q_[ff,tt]; tmp_
          m_Q_ <- tmp_
        }
      }

      for (s in 1:n) {
        params_ <- params %>% transpose() %>% pluck(s); params_
        m_Q[[s]] <- apply(m_Q_,c(1,2),function(x) eval(str2expression(x), envir=params_)); m_Q[[s]]
        m_Q[[s]] <- m_Q[[s]] * params_$n_cycle_length
        P_[[s]] <-   expm(m_Q[[s]]) %>%
          as.matrix(); P_[[s]]
      }
      P_
    })) %>%
    set_names(strategies)
  m_P %>% transpose()
}


# 2.1. Transition matrices for base case
m_P <- get_transition_matrices(strategies = v_names_str, params = params) %>% pluck(1)

