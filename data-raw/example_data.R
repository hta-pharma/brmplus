source("../compare/getProbScalarRR.R")
source("../compare/getProbScalarRD.R")
source("../compare/MyFunc.R")

get.truth.params <- function(param,
                             event = c("common", "rare"),
                             hypothesis = c("null", "alternative")) {
  event <- match.arg(event)
  hypothesis <- match.arg(hypothesis)

  if (param == "RR") {
    if (event == "common") {
      if (hypothesis == "null") {
        alpha.true <- 0
        beta.true <- c(1.5, 0.6)
        gamma.true <- c(0.2, -0.5)
      } else { # alternative
        alpha.true <- 0.3
        beta.true <- c(1.65, 0.5)
        gamma.true <- c(0.2, -0.5)
      }
    } else { # rare
      if (hypothesis == "null") {
        alpha.true <- 0
        beta.true <- c(-4.7, 0.5)
        gamma.true <- c(0.2, -0.5)
      } else { # alternative
        alpha.true <- 0.7
        beta.true <- c(-5.5, 0.5)
        gamma.true <- c(0.2, -0.5)
      }
    }
  } else if (param == "RD") {
    if (event == "common") {
      if (hypothesis == "null") {
        alpha.true <- 0
        beta.true <- c(0.9, 0.5)
        gamma.true <- c(0.2, -0.5)
      } else { # alternative
        alpha.true <- 0.1
        beta.true <- c(0.9, 0.2)
        gamma.true <- c(0.2, -0.5)
      }
    } else { # rare
      if (hypothesis == "null") {
        alpha.true <- 0
        beta.true <- c(-4.5, 0.5)
        gamma.true <- c(0.2, -0.5)
      } else { # alternative
        alpha.true <- 0.05
        beta.true <- c(-5.5, 0.2)
        gamma.true <- c(0.2, -0.5)
      }
    }
  } else {
    stop("param must be 'RR' or 'RD'")
  }

  list(
    alpha.true = alpha.true,
    beta.true = beta.true,
    gamma.true = gamma.true
  )
}


data.generate <- function(param, distribution = "unif", n, alpha.true, beta.true, gamma.true) {
  getProb <- if (param == "RR") getProbRR else getProbRD

  v.1 <- rep(1, n) # intercept term
  if (distribution == "unif") {
    v.2 <- runif(n, 0, 0.6)
  } else if (distribution == "binom") {
    v.2 <- rbinom(n, 1, 0.3)
  } else if (distribution == "norm") {
    v.2 <- rnorm(n, 0.5, 0.3)
  }
  v <- cbind(v.1, v.2)
  v.1 <- as.matrix(v.1, ncol = 1)
  pscore.true <- exp(v %*% gamma.true) / (1 + exp(v %*% gamma.true))
  p0p1.true <- getProb(v.1 %*% alpha.true, v %*% beta.true)
  x <- rbinom(n, 1, pscore.true)
  pA.true <- p0p1.true[, 1]
  pA.true[x == 1] <- p0p1.true[x == 1, 2]
  y <- rbinom(n, 1, pA.true)

  Na0 <- sum(x == 0)
  Na1 <- sum(x == 1)
  N0_1 <- sum(y[which(x == 0)])
  N1_1 <- sum(y[which(x == 1)])

  data.simulation <- list(data = data.frame(y, x, v), count = c(Na0, Na1, N0_1, N1_1))
  return(data.simulation)
}


generate.example.data <- function(
  n = 50,
  seed = 1234,
  distributions = c("unif", "binom", "norm")
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  #
  distributions <- match.arg(distributions, choices = c("unif", "binom", "norm"), several.ok = TRUE)

  # param × event × hypothesis × distribution
  cases <- expand.grid(
    param = c("RR", "RD"),
    event = c("common", "rare"),
    hypothesis = c("null", "alternative"),
    distribution = distributions,
    stringsAsFactors = FALSE
  )

  ipd_list <- vector("list", nrow(cases))
  counts_list <- vector("list", nrow(cases))

  for (i in seq_len(nrow(cases))) {
    param_i <- cases$param[i]
    event_i <- cases$event[i]
    hypothesis_i <- cases$hypothesis[i]
    distribution_i <- cases$distribution[i]

    truth <- get.truth.params(
      param      = param_i,
      event      = event_i,
      hypothesis = hypothesis_i
    )

    sim_i <- data.generate(
      param = param_i,
      distribution = distribution_i,
      n = n,
      alpha.true = truth$alpha.true,
      beta.true = truth$beta.true,
      gamma.true = truth$gamma.true
    )

    dat_i <- sim_i$data
    dat_i$param <- param_i
    dat_i$event <- event_i
    dat_i$hypothesis <- hypothesis_i
    dat_i$distribution <- distribution_i

    ipd_list[[i]] <- dat_i

    counts_list[[i]] <- data.frame(
      param = param_i,
      event = event_i,
      hypothesis = hypothesis_i,
      distribution = distribution_i,
      Na0 = sim_i$count[1],
      Na1 = sim_i$count[2],
      N0_1 = sim_i$count[3],
      N1_1 = sim_i$count[4]
    )
  }

  ipd_all <- do.call(rbind, ipd_list)
  rownames(ipd_all) <- NULL

  counts_all <- do.call(rbind, counts_list)
  rownames(counts_all) <- NULL

  list(
    ipd    = ipd_all,
    counts = counts_all
  )
}


## n = 50
exdat.50 <- generate.example.data(n = 50, seed = 1234)
example.ipd.50 <- exdat.50$ipd
example.counts.50 <- exdat.50$counts

## n = 200
exdat.200 <- generate.example.data(n = 200, seed = 1234)
example.ipd.200 <- exdat.200$ipd
example.counts.200 <- exdat.200$counts

## n = 500
exdat.500 <- generate.example.data(n = 500, seed = 1234)
example.ipd.500 <- exdat.500$ipd
example.counts.500 <- exdat.500$counts

## Output
usethis::use_data(
  example.ipd.50,
  example.counts.50,
  example.ipd.200,
  example.counts.200,
  example.ipd.500,
  example.counts.500,
  internal  = FALSE,
  overwrite = TRUE
)
