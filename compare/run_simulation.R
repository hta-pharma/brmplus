suppressPackageStartupMessages({
  library(doSNOW)
  library(doRNG)
  library(brmplus)
  library(epitools)
  library(geepack)
  library(sandwich)
  library(lmtest)
  library(brglm2)
  library(logistf)
  library(binom)
  library(epiR)
  library(PropCIs)
  library(MASS)
})

mat_vec_mul <- getFromNamespace("mat_vec_mul", "brmplus")
compute_augmentation_cpp <- getFromNamespace("compute_augmentation_cpp", "brmplus")

source("getProbScalarRR.R")
source("getProbScalarRD.R")
source("1_CallMLE.R")
source("1.1_MLE_Point.R")
source("MLE_Point_Firth_for_RR.R")
source("MLE_Point_Firth_for_RD.R")
source("1.2_MLE_Var.R")
source("bayes_p.R")
source("MyFunc.R")
source("data_generation_simulation.R")

## Defaults: one invocation covers all 24 scenarios using seeds 1:1000.
n_values <- "50,200,500"
param_values <- "RR,RD"
event_values <- "common,rare"
hypothesis_values <- "null,alternative"
result_dir <- file.path("results", "simulation")

argv <- commandArgs(TRUE)
if (length(argv) == 0L) {
  message("No arguments supplied; using the default scenario grid.")
} else {
  for (arg in argv) eval(parse(text = arg))
}

split_values <- function(x) trimws(strsplit(as.character(x), ",", fixed = TRUE)[[1]])
n_vec <- as.integer(split_values(n_values))
param_vec <- split_values(param_values)
event_vec <- split_values(event_values)
hypothesis_vec <- split_values(hypothesis_values)

if (anyNA(n_vec) || any(n_vec <= 0L)) stop("n_values must contain positive integers")
if (!all(param_vec %in% c("RR", "RD"))) stop("param_values must contain only RR and RD")
if (!all(event_vec %in% c("common", "rare"))) stop("event_values must contain only common and rare")
if (!all(hypothesis_vec %in% c("null", "alternative"))) {
  stop("hypothesis_values must contain only null and alternative")
}
seed_vec <- seq_len(1000L)
R <- length(seed_vec)

scenarios <- expand.grid(
  n = n_vec,
  event = event_vec,
  hypothesis = hypothesis_vec,
  param = param_vec,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

ncores <- min(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")), length(seed_vec), 101L)
if (is.na(ncores) || ncores < 1L) ncores <- 1L

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
cl <- makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

## foreach RNG stream registration remains disabled.
# registerDoRNG(1234)

result_by_seed <- tryCatch(
  foreach(
    seed = seed_vec,
    .packages = c(
      "brmplus", "epitools", "geepack", "sandwich", "lmtest", "brglm2",
      "MASS", "logistf", "binom", "epiR", "PropCIs"
    )
  ) %dopar% {
    lapply(seq_len(nrow(scenarios)), function(i) {
      scenario <- scenarios[i, ]
      ans <- run(
        scenario$param, scenario$n, scenario$event, scenario$hypothesis,
        seed = seed
      )
      as.data.frame(list(
        seed = seed,
        estimate = ans[1, ],
        se = ans[2, ],
        low = ans[3, ],
        up = ans[4, ],
        p = ans[5, ]
      ))
    })
  },
  finally = stopCluster(cl)
)

for (i in seq_len(nrow(scenarios))) {
  scenario <- scenarios[i, ]
  result_all <- do.call(rbind, lapply(result_by_seed, `[[`, i))
  output_file <- file.path(
    result_dir,
    paste0(
      "othermethods_results_", scenario$param, "_", scenario$event, "_",
      scenario$hypothesis, "_n_", scenario$n, "_R_", R, ".csv"
    )
  )
  write.csv(result_all, file = output_file)
}

message(
  "Completed ", length(seed_vec), " seeds x ", nrow(scenarios),
  " scenarios (1 through ", R, ")."
)
