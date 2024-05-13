###################################
# HEALTHY-SICK-DEAD
###################################
library(tictoc)
tic()
run_psa = TRUE

source(here::here("_healthy-sick-dead/01-healthy-sick-dead_setup.r"))
source(here('_healthy-sick-dead/functions_healthy-sick-dead.r'))
source(here("_healthy-sick-dead/02_healthy-sick-dead_process-parameters.r"))
if (run_psa) source(here("_healthy-sick-dead/03_healthy-sick-dead_process-parameters_psa.r"))
source(here("_healthy-sick-dead/04_healthy-sick-dead_discounting-and-cycle-adjustments.r"))
source(here("_healthy-sick-dead/05_healthy-sick-dead_construct-transition-probability-matrices.r"))
if (run_psa) source(here("_healthy-sick-dead/05b_healthy-sick-dead_construct-transition-probability-matrices-psa.r"))
source(here("_healthy-sick-dead/06_healthy-sick-dead_process-payoffs.r"))
if (run_psa) source(here("_healthy-sick-dead/06b_healthy-sick-dead_process-payoffs-psa.r"))
source(here("_healthy-sick-dead/07_healthy-sick-dead_calculate-markov-trace.r"))
if (run_psa) source(here("_healthy-sick-dead/07b_healthy-sick-dead_calculate-markov-trace-psa.r"))
source(here("_healthy-sick-dead/08_healthy-sick-dead_cea-results-base-case.r"))
if (run_psa) source(here("_healthy-sick-dead/08b_healthy-sick-dead_cea-results-psa.r"))

library(gt)

table_cea %>% gt() %>%
  sub_missing(missing_text="") %>%
  tab_source_note("ED = Dominated (Extended)\n D=Dominated")

p_cea_frontier
toc()




