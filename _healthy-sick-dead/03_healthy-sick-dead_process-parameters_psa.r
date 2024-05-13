## 1.2 Probabilistic Sensitivity Analysis

# Probability Distribution lookup table. This maps the
# distribution to the appropriate quantile function
# command in R.

dist_lut <- list("gamma" = "qgamma",
                 "lognormal" = "qlnorm",
                 "beta" = "qbeta",
                 "normal" = "qnorm",
                 "unif" = "qunif")


# 1.2.1.  Parameters with uncertainty distributions
params_psa_ <-
  params_raw  %>%
  filter(!is.na(distribution)) %>%
  select(param,distribution) %>%
  na.omit() %>%
  deframe() %>%
  as.list(); params_psa_

# 1.2.2. Halton sequence

psa_ <- halton(n=params_sc$n_sim, dim = length(params_psa_))
colnames(psa_) <- unlist(params_psa_)
colnames(psa_) <- lapply(colnames(psa_),function(x) sub_in_expression(x))

# 1.2.3. Map the Halton sequence draws to the specified distributions
# using the appropriate quantile function.

params_psa <-
  colnames(psa_) %>%
  map(~(glue("{gsub(')$','',.x)}, p={psa_[,.x]})")))  %>%
  set_names(names(params_psa_)) %>%
  bind_cols() %>% #glimpse()
  rowwise() %>%
  mutate_all(~eval(str2expression(.))) %>%
  ungroup() %>%
  as.list()

# 1.2.4. Obtain replicate copies of any parameters that do not
# have uncertainty distributions defined in the parameter
# input table.

params_sc_psa <-
  params_sc[setdiff(names(params_sc),names(params_psa))] %>%
  map(~(rep(.x,params_sc$n_sim)))

# 1.2.5. Append everything together to obtain the final
# PSA list object.

params_psa <-
  append(params_psa, params_sc_psa)

params_psa <-
  append(params_psa,map(params_fmla,~(eval(str2expression(.x),envir=params_psa))))

# 1.2.6. Ensure that the mean of the PSA draws is close to the mean of the
# base_case values.

params_psa %>%
  map(~(mean(.x))) %>%
  bind_cols() %>%
  mutate(source = "psa") %>%
  bind_rows(
    params %>% bind_cols() %>% mutate(source = "baseline")
  )  %>%
  gather(param,value,-source) %>%
  spread(source,value) %>%
  mutate(abs_diff = abs(psa-baseline)) %>%
  arrange(desc(abs_diff))

