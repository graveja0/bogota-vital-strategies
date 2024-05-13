input_file <- normalizePath(here("case-studies/healthy-sick-dead_MASTER.xlsx")); input_file
params_raw <- readxl::read_xlsx(input_file,sheet="parameters"); params_raw
payoffs_raw <- readxl::read_xlsx(input_file, sheet = "payoffs"); payoffs_raw
transition_probs_raw <- readxl::read_xlsx(input_file, sheet = "transition_probs"); transition_probs_raw
transition_rates_raw <- readxl::read_xlsx(input_file, sheet = "transition_rates"); transition_rates_raw
###########################
## Metaparameters
###########################
v_names_states <- unique(payoffs_raw$state); v_names_states
v_n_states <- length(v_names_states)
v_names_str <- unique(payoffs_raw$strategy) %>% na.omit() %>% as.vector(); v_names_str
n_strategies <- length(v_names_str)
n_sim <- params_raw %>% filter(param == "n_sim") %>% pull(base_case); n_sim

##############################
## Step 1: Process Parameters
##############################
## 1.1 BASE CASE
# 1.1.1. Scalar parameters
params_sc <-
 params_raw  %>%
 filter(is.na(formula)) %>%
 select(param,base_case) %>%
 na.omit() %>%
 deframe() %>%
 as.list(); params_sc
# 1.1.2. Parameters that are functions of other parameters
params_fmla <-
 params_raw %>%
 filter(!is.na(formula)) %>%
 select(param,formula) %>%
 na.omit() %>%
 deframe() %>%
 as.list(); params_fmla
# 1.1.3. Execute the formulas and create the master params object
params <-
 append(params_sc,map(params_fmla,~(eval(str2expression(.x),envir=params_sc))))
# Adding another meta-parameter that is a function.
n_cycles <- params$n_cycles
