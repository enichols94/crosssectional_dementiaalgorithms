##########################################################################
### Author: Emma Nichols
### Date: 1/26/2022
### Purpose: KM Curves for Test set
##########################################################################

rm(list = ls())

pacman::p_load(data.table, openxlsx, haven, readr, dplyr, magrittr, stringr,
               labelled, scales, gridExtra, grid, ggplot2, mice, survival,
               survminer, arules, gglasso, cowplot, pROC, survival, survminer, 
               cutpointr)
date <- gsub("-", "_", Sys.Date())

cs_model_date <- DATE

# SET OBJECTS -------------------------------------------------------------

dir <- FILEPATH
rawdata_dir <- paste0(dir, "raw_data/")
derived_dir <- paste0(dir, "derived_data/")
plot_dir <- paste0(dir, "paper/km_plot/")

# GET DATA ----------------------------------------------------------------

variable_map <- as.data.table(read.xlsx(paste0(dir, "variable_map.xlsx")))
alg_items <- variable_map[, new_name]

# GET PREDICTIONS ---------------------------------------------------------

get_preddata <- function(date, model = "model"){
  model_dir <- paste0(FILEPATH, date, "/")
  model <- read_rds(paste0(model_dir, model, ".rds"))
  
  dt <- readr::read_rds(paste0(derived_dir, "rosmap_formatted.rds"))
  
  model_names <- gsub("[0-9]$", "", rownames(model$gglasso.fit$beta))
  alg_items <- alg_items[grepl(paste0(model_names, collapse = "|"), alg_items)]
  
  # FORMAT VARIABLES --------------------------------------------------------
  
  message("Formatting data")
  
  ## BINARY - CREATE MISSING INDICATOR
  ## <15 CATEGORIES - CREATE MODEL MATRIX AND MISSING INDICATOR
  ## >=15 CATEGORIES - EQUAL INTERVAL DISCRETIZATION AND MISSING INDICATOR
  
  for (var in c(alg_items)){
    nunique <- dt[!is.na(get(var)), length(unique(get(var)))]
    if (nunique == 2){ ## binary
      dt[, paste0(var, "_ind") := as.numeric(is.na(get(var)))]
      dt[is.na(get(var)), c(var) := 0]
    } else if (nunique > 2 & nunique < 15){ ## factor no collapse
      dt[, c(var) := as.factor(get(var))]
      
      ## relevel so largest factor is reference
      levels <- dt[!is.na(get(var)), .N, by = var]
      ref_level <- as.character(levels[N == max(N), get(var)])
      dt[, c(var) := relevel(get(var), ref = ref_level)]
      
      ## get model matrix with NA's 
      mat <- as.data.table(model.matrix(~., data = model.frame(~., as.data.frame(dt[, c(var), with = F]), na.action=na.pass))[, -1]) 
      mat[, c(names(mat)) := lapply(.SD, function(x) ifelse(is.na(x), 0, x)), .SDcols = names(mat)] ## set NA's to 0
      mat[, paste0(var, "_ind") := dt[, as.numeric(is.na(get(var)))]]
      mat[, projid := dt[, projid]]; mat[, fu_year := dt[, fu_year]] ## attach on ids
      dt <- merge(dt, mat, by = c("projid", "fu_year"), sort = F)
      dt[, c(var) := NULL]
    } else if (nunique >=15){ ## factor with collapse
      dt[, c(var) := discretize(get(var), breaks = 14, method = "interval", labels = F)]
      dt[, c(var) := as.factor(get(var))]
      
      ## relevel so largest factor is reference
      levels <- dt[!is.na(get(var)), .N, by = var]
      ref_level <- as.character(levels[N == max(N), get(var)])
      dt[, c(var) := relevel(get(var), ref = ref_level)]
      
      ## get model matrix with NA's 
      mat <- as.data.table(model.matrix(~., data = model.frame(~., as.data.frame(dt[, c(var), with = F]), na.action=na.pass))[, -1]) 
      mat[, c(names(mat)) := lapply(.SD, function(x) ifelse(is.na(x), 0, x)), .SDcols = names(mat)] ## set NA's to 0
      mat[, paste0(var, "_ind") := dt[, as.numeric(is.na(get(var)))]]
      mat[, projid := dt[, projid]]; mat[, fu_year := dt[, fu_year]] ## attach on ids
      dt <- merge(dt, mat, by = c("projid", "fu_year"), sort = F)
      dt[, c(var) := NULL]
    }
  }
  
  # GET OTHER VARS READY ----------------------------------------------------
  
  dt[, `:=` (raceother = as.numeric(race == "other"), raceblack = as.numeric(race == "black"))]
  
  ## EDUC
  dt[, educ_u12 := as.numeric(educ < 12)][, educ_o1216 := as.numeric(educ > 12 & educ <= 16)]
  dt[, educ_1718 := as.numeric(educ >= 17 & educ <= 18)][, educ_19 := as.numeric(educ > 18)]
  
  ## AGE
  dt[, age := discretize(age_at_visit, breaks = 14, method = "interval", labels = F)]
  dt[age %in% c(1,2,4), age := 5][age == 14, age := 13] ## take in the tails a bit
  dt[, age := as.factor(age)]
  ## relevel so largest factor is reference
  levels <- dt[!is.na(age), .N, by = "age"]
  ref_level <- as.character(levels[N == max(N), age])
  dt[, age := relevel(age, ref = ref_level)]
  ## get model matrix with NA's 
  mat <- as.data.table(model.matrix(~., data = model.frame(~., as.data.frame(dt[, "age", with = F]), na.action=na.pass))[, -1]) 
  age_names <- names(mat)[grepl("age", names(mat))] ## get names of binary age vars
  mat[, projid := dt[, projid]]; mat[, fu_year := dt[, fu_year]] ## attach on ids
  dt <- merge(dt, mat, by = c("projid", "fu_year"), sort = F)
  
  otheralg_items <- c(age_names, "msex", "educ_u12", "educ_o1216", "educ_1718", "educ_19", "raceblack", "raceother", "hispanic")
  
  # GET TEST DATA -----------------------------------------------------------
  
  test_dt <- copy(dt[testdata == 1])
  
  message("Getting predictions")
  
  alg_preds <- sort(names(test_dt)[grepl(paste0(alg_items, collapse = "|"), names(test_dt))])
  alg_preds <- c(otheralg_items, alg_preds)

  test_dt[, pred := as.numeric(predict(model, newx = as.matrix(test_dt[, c(alg_preds), with = F]), s = "lambda.min") == 1)]

  
  return(test_dt)
}

fullmodel_preds <- get_preddata(cs_model_date)
smallmodel_preds <- get_preddata(cs_model_date, model = "smallest_model")

# GET SURVIVAL DATA -------------------------------------------------------

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
  
  ## subset to non-missing records
  subset_dt <- copy(dt[, c("projid", "fu_year", outcome, "age_at_visit", "fu_max_i", "diabetes_tv", "age_bl"), with = F])
  subset_dt <- na.omit(subset_dt)
  
  ## number excluded
  print(paste0("Number excluded: ", dt[, length(unique(projid))] - subset_dt[, length(unique(projid))]))
  
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

sdt_goldstandard <- do_survival(fullmodel_preds, "dementia")
sdt_full <- do_survival(fullmodel_preds, "pred")
sdt_small <- do_survival(smallmodel_preds, "pred")

# CREATE GRAPHS -----------------------------------------------------------

# PANEL A - GOLD STANDARD -------------------------------------------------

sfit_gs <- survfit(Surv(time, end_time, outcome) ~ diabetes_tv, data = sdt_goldstandard)
splot_gs <- ggsurvplot(sfit_gs, palette = c("deepskyblue4", "deepskyblue3"), 
                       legend.title = "", legend.labs = c("No Diabetes", "Diabetes"),
                       xlab = "Time (years from baseline)", title = "A.", legend = "bottom")

# PANEL B FULL ALGORITHM PREV/INC -----------------------------------------

sfit_f <- survfit(Surv(time, end_time, outcome) ~ diabetes_tv, data = sdt_full)
sfit_s <- survfit(Surv(time, end_time, outcome) ~ diabetes_tv, data = sdt_small)

splot_fcomb <- ggsurvplot(fit = list("Gold Standard" = sfit_gs, "Full Algorithm" = sfit_f, "Small Algorithm" = sfit_s), 
                          combine = T, legend.title = "", 
                          legend.labs = c("Gold Standard:\nNo Diabetes", "Gold Standard:\nDiabetes",
                                          "Full Algorithm:\nNo Diabetes", "Full Algorithm:\nDiabetes",
                                          "Small Algorithm:\nNo Diabetes", "Small Algorithm:\nDiabetes"),
                          xlab = "Time (years from baseline)", title = "B.", legend = "bottom",
                          palette = c("deepskyblue4", "deepskyblue3", 
                                      "mediumpurple4", "mediumpurple3",
                                      "#F29D7E", "#F22E2E")) +
  guides(colour = guide_legend(nrow = 2))

full_graph <- plot_grid(splot_gs$plot, splot_fcomb$plot, nrow = 1, align = "h", axis = "tb")

ggsave(paste0(plot_dir, "km_plot_", date, ".tiff"), plot = full_graph, height = 5, width = 10)
