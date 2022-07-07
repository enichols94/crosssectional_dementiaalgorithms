##########################################################################
### Author: Emma Nichols
### Date: 1/26/2022
### Purpose: Simulation study
##########################################################################

rm(list = ls())

pacman::p_load(data.table, openxlsx, haven, readr, dplyr, magrittr, stringr,
               labelled, scales, gridExtra, grid, ggplot2, mice, survival,
               survminer, arules, gglasso, cowplot, pROC, survival, survminer,
               simsurv, flexsurv, pbapply, coxed)
date <- gsub("-", "_", Sys.Date())
set.seed(12494)

# GET ARGUMENTS -----------------------------------------------------------

## SETUP TO PASS CENSOR TIME AS AN ARGUMENT TO SYSTEM
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
censor_time <- as.numeric(args[1]) ## censor time - either 3, 5, 10, 20
print(censor_time)

# SET OBJECTS -------------------------------------------------------------

dir <- FILEPATH
rawdata_dir <- paste0(dir, "raw_data/")
derived_dir <- paste0(dir, "derived_data/")
plot_dir <- paste0(dir, "plots/")

# GET DATA ----------------------------------------------------------------

rosmap_dt <- readr::read_rds(paste0(derived_dir, "rosmap_formatted.rds"))

# FORMAT SURVIVAL DATA ----------------------------------------------------

do_survival <- function(data, outcome){
  
  dt <- copy(data)
  dt <- dt[!is.na(get(outcome))]
  
  ## get rid of prevalent outcome
  dt <- dt[order(projid, fu_year)]
  suppressWarnings(dt[, any_dementia := as.numeric(max(get(outcome), na.rm = T)), by = "projid"])
  dt[is.infinite(any_dementia), any_dementia := 0]
  dt[get(outcome) == 1, first_dem_visit := min(fu_year), by = "projid"]
  dt[any_dementia == 1, first_dem_visit := max(first_dem_visit, na.rm = T), by = "projid"] ## set first dementia visit for all visits among those with dementia
  dt[is.na(first_dem_visit), first_dem_visit := 100] ## set to arbitrarily large and non-existent visit
  dt <- dt[fu_year <= first_dem_visit]
  
  ## define max visit 
  dt[!is.na(diabetes_tv), fu_max := max(fu_year), by = "projid"]
  dt[, fu_max_i := as.numeric(fu_year == fu_max)]
  
  ## indicator for if you ever got dementia
  dt[, case_i := max(get(outcome)), by = "projid"]
  
  # MISSING DATA PATTERNS ---------------------------------------------------
  
  covs <- c("smoke_cf", "hyp_bl", "claud_bl", "heart_bl", "stroke_bl", "bmi_bl", "diabetes_tv",
            "race", "hispanic", "study", "educ", "msex", "age_bl")
  
  ## subset to non-missing records
  subset_dt <- copy(dt[, c("projid", "fu_year", outcome, "age_at_visit", "fu_max_i", covs), with = F])
  subset_dt <- na.omit(subset_dt)
  
  ## reset follow-up year to be in units
  subset_dt <- subset_dt[order(projid, fu_year)]
  subset_dt[, fu_year := sequence(.N)-1, by = "projid"]
  
  # SET UP SURVIVAL ANALYSIS ------------------------------------------------
  
  ## get time scale
  subset_dt[, time := (age_at_visit - age_bl)]
  
  ## get table of incident dementia status
  get_incdem <- function(wave){
    inc_wave <- wave+1
    inc_dementia <- copy(subset_dt[fu_year == inc_wave, .(projid, outcome = get(outcome), end_time = time, 
                                                          fu_year = wave)])
    return(inc_dementia)
  }
  inc_dem_dt <- rbindlist(lapply(0:(dt[, max(fu_year)]-1), get_incdem))
  
  ## merge back on subset 
  subset_dt[, c(outcome) := NULL]
  subset_dt <- merge(subset_dt, inc_dem_dt, all.x = T, by = c("projid", "fu_year"))
  
  ## drop last visit
  subset_dt <- subset_dt[!fu_max_i == 1]
  
  return(subset_dt)
}

rosmap_survdt <- do_survival(rosmap_dt, "dementia")

# GET INPUTS --------------------------------------------------------------

flexmod <- flexsurv::flexsurvspline(Surv(end_time, outcome) ~ 1, data = rosmap_survdt, k = 3)

plot(flexmod, 
     main = "Flexible parametric model",
     ylab = "Survival probability",
     xlab = "Time")

cox_reg <- coxph(Surv(time, end_time, outcome) ~ age_bl + I(age_bl^2) + msex + educ + race + hispanic + as.factor(study) + diabetes_tv + 
                   stroke_bl + claud_bl + smoke_cf + heart_bl + bmi_bl + hyp_bl, data = rosmap_survdt)

## get sample size
rosmap_survdt[!is.na(age_bl) & !is.na(msex) &!is.na(educ) & !is.na(race) & !is.na(hispanic) &!is.na(study) & !is.na(diabetes_tv) & !is.na(stroke_bl) & 
                !is.na(claud_bl) & !is.na(smoke_cf) & !is.na(heart_bl) & !is.na(bmi_bl) & !is.na(hyp_bl), length(unique(projid))]

rosmap_dt[fu_year == 0, mean(diabetes_tv, na.rm = T)]

# SIMULATE DATA -----------------------------------------------------------

# Define a function returning the log cum hazard at time t - will need to alter
logcumhaz <- function(t, x, betas, knots) {
  
  # Obtain the basis terms for the spline-based log
  # cumulative hazard (evaluated at time t)
  basis <- flexsurv::basis(knots, log(t))
  
  # Evaluate the log cumulative hazard under the
  # Royston and Parmar specification
  res <- 0
  for (n in 1:length(knots)){
    res <- res + betas[[paste0("gamma", n-1)]] * basis[[n]]
  }
  res <- res + betas[["diabetes"]] * x[["diabetes"]]
  
  # Return the log cumulative hazard at time t
  res
}

get_performance <- function(row, expanded_dt, c, iter){
  
  message(paste0(iter, " - Stat: ", row))
  
  sens <- smallparam_dt[, sens][row]
  spec <- smallparam_dt[, spec][row]
  ratio_dsens <- smallparam_dt[, ratio_dsens][row]
  ratio_dspec <- smallparam_dt[, ratio_dspec][row]
  
  alt_dt <- copy(expanded_dt)
  alt_dt[, max_fu := max(fu_year), by = "id"]
  
  ## GENERATE RANDOM ALGORITHM PERFORMANCE BY SENS/SPEC
  alt_dt[status == 0 & diabetes == 0, alt_status := rbinom(nrow(alt_dt[status == 0 & diabetes == 0]),1,(1-spec))]
  alt_dt[status == 1 & diabetes == 0, alt_status := rbinom(nrow(alt_dt[status == 1 & diabetes == 0]),1,(sens))]
  alt_dt[status == 0 & diabetes == 1, alt_status := rbinom(nrow(alt_dt[status == 0 & diabetes == 1]),1,(1-spec*ratio_dspec))]
  alt_dt[status == 1 & diabetes == 1, alt_status := rbinom(nrow(alt_dt[status == 1 & diabetes == 1]),1,(sens)*ratio_dsens)]
  ## GO BACK AND CONTINUE FOR THOSE THAT WE HAVEN'T "CAPTURED YET" - AS IF CONTINUED
  ## TO BE OBSERVED
  alt_dt[, pred_case := max(alt_status), by = "id"]
  alt_dt[, true_case := max(status), by = "id"]
  alt_dt[, unfinished := as.numeric(!max_fu == c & !pred_case == true_case)] ## unfinished if pred case is not equal true case and haven't reached max follow-up
  if (nrow(alt_dt[unfinished == 1]) > 0){
    alt_expand_dt <- copy(alt_dt[unfinished == 1])
    alt_expand_dt[, nexpand := c - max_fu]
    alt_expanded_dt <- alt_expand_dt[fu_year==max_fu][rep(1:.N,nexpand)][,fu_year:=1:.N,by=id]
    alt_expanded_dt[, fu_year := fu_year + max_fu]
    alt_expanded_dt[status == 0 & diabetes == 0, alt_status := rbinom(nrow(alt_expanded_dt[status == 0 & diabetes == 0]),1,(1-spec))]
    alt_expanded_dt[status == 1 & diabetes == 0, alt_status := rbinom(nrow(alt_expanded_dt[status == 1 & diabetes == 0]),1,(sens))]
    alt_expanded_dt[status == 0 & diabetes == 1, alt_status := rbinom(nrow(alt_expanded_dt[status == 0 & diabetes == 1]),1,(1-spec*ratio_dspec))]
    alt_expanded_dt[status == 1 & diabetes == 1, alt_status := rbinom(nrow(alt_expanded_dt[status == 1 & diabetes == 1]),1,(sens)*ratio_dsens)]
    final_alt_dt <- rbindlist(list(alt_dt[unfinished == 0], alt_expand_dt, alt_expanded_dt),
                              fill = T, use.names = T)
    rm(list = c("alt_dt", "alt_expand_dt", "alt_expanded_dt"))
  } else {
    final_alt_dt <- copy(alt_dt)
    rm(list = c("alt_dt"))
  }
  ## FIRST DEMENTIA VISIT IS INCIDENT 
  final_alt_dt[alt_status == 1, pred_fv := min(fu_year), by = "id"]
  suppressWarnings(final_alt_dt[, pred_fv := max(pred_fv, na.rm = T), by = "id"])
  final_alt_dt[is.na(pred_fv), pred_fv := 100] ## arbitrarily large
  final_alt_dt <- final_alt_dt[!fu_year > pred_fv]
  to_delete <- intersect(c("replicates", "max_fu", "pred_case", "true_case", "unfinished", "nexpand", "pred_fv"), names(final_alt_dt))
  final_alt_dt[, c(to_delete) := NULL]
  final_alt_dt[, max_fu := max(fu_year), by = "id"]
  final_alt_dt <- final_alt_dt[fu_year == max_fu]
  
  # RUN MODELS --------------------------------------------------------------
  
  ## TRUE MODEL
  true_dt <- copy(expanded_dt)
  true_dt[, fu_year_maxi := as.numeric(fu_year == max(fu_year)), by = c("id")]
  true_dt <- true_dt[fu_year_maxi == 1]
  true_reg <- coxph(Surv(fu_year, status) ~ diabetes, data = true_dt)
  
  ## ALT MODEL
  alt_reg <- coxph(Surv(fu_year, alt_status) ~ diabetes, data = final_alt_dt)
  
  rm(list = c("true_dt", "final_alt_dt"))
  
  return_dt <- data.table(sensitivity = sens, specificity = spec, sens_ratio = ratio_dsens, spec_ratio = ratio_dspec, 
                          true = as.numeric(coef(true_reg)), alt = as.numeric(coef(alt_reg)))
  return(return_dt)
  
}


simulation <- function(n, flexible_model, cov_beta, diabetes_proportion, hazard_fun = logcumhaz,
                       ct = 15, param_dt){
  print(paste0("Iteration ", n))
  
  # SIMULATE DATA -----------------------------------------------------------

  cov_dt <- data.table(id = 1:4000, diabetes = rbinom(4000, 1, diabetes_proportion))
  
  ## add in coefficient into flexmod coefficients
  flexible_model$coefficients['diabetes'] <- cov_beta
  
  sim_dt <- simsurv(betas = flexible_model$coefficients, # "true" parameter values
                 x = cov_dt,                          # covariate data
                 knots = flexible_model$knots,        # knot locations for splines
                 logcumhazard = hazard_fun,           # definition of log cum hazard
                 maxt = ct,                           # no right-censoring
                 )           
  
  sim_dt <- merge(as.data.table(sim_dt), cov_dt, by = "id")
  
  ## BUILD OUT LONGITUDINAL DATA SO CAN VARY ALGORITHM
  sim_dt[, replicates := ceiling(eventtime)]
  edt <- sim_dt[rep(1:.N,replicates)][,fu_year:=1:.N,by=id]
  edt[status == 1 & !fu_year == ceiling(eventtime), status := 0]
  rm(sim_dt)
  
  # GET RESULTS -------------------------------------------------------------

  result_dt <- rbindlist(lapply(1:nrow(param_dt), function(x) get_performance(row = x, expanded_dt = edt, c = ct, iter = n))) 
  return(result_dt)
  
 }

smallparam_dt <- as.data.table(expand.grid(sens = seq(0.4,1,0.1), spec = c(0.8,0.85,0.9,0.95,0.975,1),
                                           ratio_dsens = c(0.85,0.9,0.95,1,1.05,1.1,1.15), ratio_dspec = c(0.85,0.9,0.95,1,1.05,1.1,1.15)))
smallparam_dt <- smallparam_dt[!(spec*ratio_dspec>1) & !(sens*ratio_dsens > 1)]

print(paste0("Results for censoring time: ", censor_time))
result <- lapply(1:500, function(x) simulation(n = x, flexible_model = flexmod, cov_beta = as.numeric(coef(cox_reg)['diabetes_tv']),
                                               diabetes_proportion = 0.14, hazard_fun = logcumhaz, ct = censor_time,
                                               param_dt = smallparam_dt))

readr::write_rds(result, paste0(derived_dir, "simulations_", censor_time, "_", date, ".rds"))
readr::write_rds(smallparam_dt, paste0(derived_dir, "parameters_", censor_time, "_", date, ".rds"))

