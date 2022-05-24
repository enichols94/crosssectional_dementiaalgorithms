##########################################################################
### Author: Emma Nichols
### Date: 8/12/2021
### Purpose: Explore baseline survival analysis
##########################################################################

rm(list = ls())

pacman::p_load(data.table, openxlsx, haven, readr, dplyr, magrittr, stringr,
               labelled, scales, gridExtra, grid, ggplot2, mice, survival,
               survminer)
date <- gsub("-", "_", Sys.Date())

# SET OBJECTS -------------------------------------------------------------

dir <- FILEPATH
rawdata_dir <- paste0(dir, "raw_data/")
derived_dir <- paste0(dir, "derived_data/")

# GET DATA ----------------------------------------------------------------

dt <- readr::read_rds(paste0(derived_dir, "rosmap_formatted.rds"))

## get rid of prevalent dementia
dt <- dt[prev_dementia == 0 | is.na(prev_dementia)]

## define max visit
dt[!is.na(diabetes_tv), fu_max := max(fu_year), by = "projid"] 

## indicator for if you ever got dementia
dt[, case_i := max(inc_dementia), by = "projid"]

## indicator for missing dementia status
dt[, missing_dem := as.numeric(is.na(dementia))]

# MISSING DATA PATTERNS ---------------------------------------------------

covs <- c("smoke_cf", "hyp_bl", "claud_bl", "heart_bl", "stroke_bl", "bmi_bl", "diabetes_tv",
          "race", "hispanic", "study", "educ", "msex", "age_bl")

md.pattern(dt[, c(covs), with = F], plot = F)

## subset to non-missing records
subset_dt <- copy(dt[, c("projid", "fu_year", "inc_dementia", "age_at_visit", "fu_max_i",
                         "missing_dem", covs), with = F])
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
  inc_dementia <- copy(subset_dt[fu_year == inc_wave, .(projid, dementia = inc_dementia, end_time = time, 
                                                        fu_year = wave, missing_dem)])
  return(inc_dementia)
}
inc_dem_dt <- rbindlist(lapply(0:(dt[, max(fu_year)]-1), get_incdem))

## merge back on subset 
subset_dt[, c("missing_dem", "inc_dementia") := NULL]
subset_dt <- merge(subset_dt, inc_dem_dt, all.x = T, by = c("projid", "fu_year"))

## drop last visit
subset_dt <- subset_dt[!fu_max_i == 1]

## drop when dementia is missing at the last visit
subset_dt[, last_interval := max(fu_year), by = "projid"]
subset_dt <- subset_dt[(last_interval == fu_year & missing_dem == 0) | !last_interval == fu_year] 

# DO SURVIAL ANALYSIS -----------------------------------------------------

## center variables
subset_dt[, agec_bl := age_bl - 75]

## survival curve
sfit <- survfit(Surv(time, end_time, dementia) ~ diabetes_tv, data = subset_dt)
ggsurvplot(sfit) ## unadjusted

## adjusted for demographics
model1 <- coxph(Surv(subset_dt[, time], subset_dt[, end_time], subset_dt[, dementia]) ~
                 agec_bl + msex + educ + race + hispanic + as.factor(study) + diabetes_tv, data = subset_dt)

## full adjusted
model_final <- coxph(Surv(time, end_time, dementia) ~
                 agec_bl + I(agec_bl^2) + msex + educ + race + hispanic + as.factor(study) + diabetes_tv + 
                 stroke_bl + claud_bl + smoke_cf + heart_bl + bmi_bl + hyp_bl, data = subset_dt)


unadjusted_model <- coxph(Surv(time, end_time, dementia) ~ diabetes_tv, data = subset_dt)

## Schoenfeld residuals
resids <- cox.zph(model_final)

plot(resids, var = "bmi_bl")

cox.zph(model_final)

## New dataset
subset2_dt <- survSplit(Surv(time, end_time, dementia) ~ ., data= subset_dt, cut=c(5,12),
                       episode= "tgroup", id="id")

model_timevar <- coxph(Surv(time, end_time, dementia) ~
                       agec_bl + I(agec_bl^2) + msex + educ + race + hispanic + as.factor(study) + 
                         diabetes_tv + stroke_bl + claud_bl + smoke_cf + heart_bl + bmi_bl:strata(tgroup) + 
                         hyp_bl:strata(tgroup), data = subset2_dt)

cox.zph(model_timevar)
