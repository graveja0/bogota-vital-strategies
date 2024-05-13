# 3.0 Function to obtain payoff vectors given specified payoffs by strategy
# as specified in Excel document.

get_payoffs <- function(params, jumpover = NULL) {
  n <- params %>% map_dbl(~(length(.x))) %>% max()

  #tic()
  out <-
    map2_progress(rep(list(payoffs_raw),n),transpose(params),~({
      tmp_ <- .y

      payoff_ <-
        .x %>%
        rowwise() %>%
        mutate_at(vars(contains("fmla")),function(x) eval(str2expression(x),envir=.y)) %>%
        rename_at(vars(contains("fmla")),function(x) gsub("_fmla","",x)) %>%
        select(strategy,state,costs,qalys)

      if (!is.null(jumpover)) {
        for (jj in names(jumpover)) {
          ff <- jumpover[[jj]][[1]]; ff
          tt <- jumpover[[jj]][[2]]; tt

          payoff_j <- payoff_ %>% filter(label == tt) %>%
            mutate(label = jj,
                   description = jj)

          payoff_ <- payoff_ %>% bind_rows(payoff_j)

          v_names_states_ <- c(v_names_states,jj)
        }
      } else {
        v_names_states_ <- v_names_states
      }

      # payoff_ %>%
      #   group_by(strategy) %>%
      #   nest() %>%
      #   deframe() %>%
      #   map(~(
      #     .x %>%
      #       select(state,costs,qalys)  %>%
      #       gather(type,value,-state) %>%
      #       group_by(type) %>%
      #       nest() %>%
      #       deframe() %>%
      #       map(~(.x %>%
      #               spread(state,value) %>%
      #               select_at(vars(v_names_states_)) %>%
      #               as.matrix()))
      #   ))

#
      payoff_ %>% gather(outcome,value,-strategy,-state) %>% spread(state,value) %>% group_by(strategy) %>%
        split(.,.$strategy)   %>%
        map(~(.x  %>% split(.,.$outcome))) %>%
        map(~(.x %>% map(~(.x %>% ungroup() %>% select(-strategy,-outcome) %>% as.matrix() %>% {.[,v_names_states] %>% as.matrix() %>% t()}))))


    }))
  #toc()
  if (n==1) out <- out %>% pluck(1)
  return(out)
}


# 3.1. Payoffs for base case.
payoffs <- get_payoffs(params)

