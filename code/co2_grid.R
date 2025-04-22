library(BayesGP)
library(tidyverse)
library(npreg)
library(parallel)

function_path <- "./code"
output_path <- "./output/co2"
data_path <- "./data/co2"
source(paste0(function_path, "/00_BOSS.R"))


### Read in the full data:
co2s = read.table(paste0(data_path, "/daily_flask_co2_mlo.csv"), header = FALSE, sep = ",",
                  skip = 69, stringsAsFactors = FALSE, col.names = c("day",
                                                                     "time", "junk1", "junk2", "Nflasks", "quality",
                                                                     "co2"))
co2s$date = strptime(paste(co2s$day, co2s$time), format = "%Y-%m-%d %H:%M",
                     tz = "UTC")
co2s$day = as.Date(co2s$date)
timeOrigin = as.Date("1960/03/30")
co2s$timeYears = round(as.numeric(co2s$day - timeOrigin)/365.25,
                       3)
co2s$dayInt = as.integer(co2s$day)
allDays = seq(from = min(co2s$day), to = max(co2s$day),
              by = "7 day")
observed_dataset <- co2s %>% filter(!is.na(co2s$co2)) %>% dplyr::select(c("co2", "timeYears"))
observed_dataset$quality <- ifelse(co2s$quality > 0, 1, 0)
observed_dataset <- observed_dataset %>% filter(quality == 0)
### Create the covariate for trend and the 1-year seasonality
observed_dataset$t1 <- observed_dataset$timeYears
observed_dataset$t2 <- observed_dataset$timeYears
observed_dataset$t3 <- observed_dataset$timeYears
observed_dataset$t4 <- observed_dataset$timeYears


lower <- 2; upper <- 6; noise_var <- 1e-6
eval_once <- function(a){
  mod_once <- model_fit(formula = co2 ~ f(x = t1, model = "IWP", order = 2, sd.prior = list(param = 30, h = 10), k = 30, initial_location = median(observed_dataset$t1)) +
                          f(x = t2, model = "sGP", sd.prior = list(param = 1, h = 10), a = (2*pi), m = 2, k = 30) +
                          f(x = t3, model = "sGP", sd.prior = list(param = 1, h = 10), a = (2*pi/a), m = 1, k = 30),
                        data = observed_dataset,
                        control.family = list(sd.prior = list(param = 1)),
                        family = "Gaussian", aghq_k = 3)
  mod_once$mod$normalized_posterior$lognormconst
}



x_vals  <- seq(2.01, 6, by = 0.01)
n_cores <- 12
# this returns a list of length length(x_vals)
res_list  <- mclapply(x_vals, eval_once, mc.cores = n_cores)
exact_vals <- unlist(res_list)

# Close the progress bar
exact_grid_result <- data.frame(x = x_vals, exact_vals = exact_vals)
exact_grid_result$exact_vals <- exact_grid_result$exact_vals - max(exact_grid_result$exact_vals)
exact_grid_result$fx <- exp(exact_grid_result$exact_vals)

# Calculate the differences between adjacent x values
dx <- diff(exact_grid_result$x)
# Compute the trapezoidal areas and sum them up
integral_approx <- sum(0.5 * (exact_grid_result$fx[-1] + exact_grid_result$fx[-length(exact_grid_result$fx)]) * dx)
exact_grid_result$pos <- exact_grid_result$fx / integral_approx
save(exact_grid_result, file = paste0(output_path, "/exact_grid_result.rda"))
