##########################################################################
### Author: Emma Nichols
### Date: 1/20/2022
### Purpose: Algorithm performance
##########################################################################

rm(list = ls())

pacman::p_load(data.table, openxlsx, haven, readr, dplyr, magrittr, stringr,
               labelled, scales, gridExtra, grid, ggplot2, mice, survival,
               survminer, arules, gglasso, cowplot, pROC, survival, survminer, 
               epiR, cutpointr, ResourceSelection)
date <- gsub("-", "_", Sys.Date())

model_date <- DATE 

# SET OBJECTS -------------------------------------------------------------

dir <- FILEPATH
rawdata_dir <- paste0(dir, "raw_data/")
derived_dir <- paste0(dir, "derived_data/")
model_dir <- paste0(FILEPATH, model_date, "/")
plot_dir <- paste0(dir, "plots/")

# GET MODELS --------------------------------------------------------------

model <- read_rds(paste0(model_dir, "model.rds"))
model2 <- read_rds(paste0(model_dir, "smaller_model.rds"))
model3 <- read_rds(paste0(model_dir, "smallest_model.rds"))

# GET DATA ----------------------------------------------------------------

dt <- readr::read_rds(paste0(derived_dir, "rosmap_formatted.rds"))

variable_map <- as.data.table(read.xlsx(paste0(dir, "variable_map.xlsx")))
alg_items <- variable_map[, new_name]
model_names <- gsub("[0-9]$", "", rownames(model$gglasso.fit$beta))
alg_items <- alg_items[grepl(paste0(model_names, collapse = "|"), alg_items)]

# FORMAT VARIABLES --------------------------------------------------------

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

## GET AGE GROUP VARIABLE
dt[, age_o75 := as.numeric(age_bl > 75)]

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

# GET TEST DT AND PREDICTIONS ---------------------------------------------

test_dt <- copy(dt[testdata == 1])

## look at sample size for various follow-ups (prevalence test sample)
test_dt[, .N, by = "fu_year"]

## all incident and no-dementia cases
test_dt_inc <- copy(test_dt[prev_dementia == 0])

## get sample size
test_dt_inc[, length(unique(projid))]

get_performance <- function(data, model){
  dt <- copy(data)
  model_names <- gsub("[0-9]$", "", rownames(model$gglasso.fit$beta))
  alg_items <- alg_items[grepl(paste0(model_names, collapse = "|"), alg_items)]
  alg_preds <- sort(names(test_dt)[grepl(paste0(alg_items, collapse = "|"), names(test_dt))])
  alg_preds <- c(otheralg_items, alg_preds)
  
  dt[, pred := as.numeric(predict(model, newx = as.matrix(dt[, c(alg_preds), with = F]), s = "lambda.min") == 1)]
  dt[, pred_logodds := predict(model, newx = as.matrix(dt[, c(alg_preds), with = F]), s = "lambda.min", type = "link")]
  dt[, pred_prob := exp(pred_logodds)/(1+exp(pred_logodds))]
  sens <- nrow(dt[dementia == 1 & pred == 1])/nrow(dt[dementia == 1])
  spec <- nrow(dt[dementia == 0 & pred == 0])/nrow(dt[dementia ==0])
  hl <- hoslem.test(dt[, dementia], dt[, pred_prob])
  roc_result <- pROC::roc(dt[, dementia], dt[, pred_prob])
  return(list = c(sensitivity = sens, specificity = spec, roc = roc_result, hl_test = hl))
}

get_roccuts <- function(data, model){
  dt <- copy(data)
  model_names <- gsub("[0-9]$", "", rownames(model$gglasso.fit$beta))
  alg_items <- alg_items[grepl(paste0(model_names, collapse = "|"), alg_items)]
  alg_preds <- sort(names(test_dt)[grepl(paste0(alg_items, collapse = "|"), names(test_dt))])
  alg_preds <- c(otheralg_items, alg_preds)
  
  dt[, pred := as.numeric(predict(model, newx = as.matrix(dt[, c(alg_preds), with = F]), s = "lambda.min") == 1)]
  dt[, pred_logodds := predict(model, newx = as.matrix(dt[, c(alg_preds), with = F]), s = "lambda.min", type = "link")]
  dt[, pred_prob := exp(pred_logodds)/(1+exp(pred_logodds))]
  
  get_cut <- function(cut){
    preds <- dt[, pred_cut := as.numeric(pred_prob > cut)]
    result_row <- data.table(cut = cut, sens = nrow(dt[dementia == 1 & pred_cut == 1])/nrow(dt[dementia == 1]),
                             spec_m1 = 1-nrow(dt[dementia == 0 & pred_cut == 0])/nrow(dt[dementia ==0]))
  }
  result_dt <- rbindlist(lapply(seq(0,1,0.001), get_cut))
  return(result_dt)
}

# TABLE OF PERFORMANCE ----------------------------------------------------

get_chunk <- function(model, name){
  results <- data.table()
  for (wave in 5:10){
    ## follow-up year
    prev_dt <- copy(test_dt[fu_year == wave]) 
    ## all incident and no-dementia cases
    inc_dt <- copy(test_dt[prev_dementia == 0])
    results_prev <- get_performance(prev_dt, model)
    results_inc <- get_performance(inc_dt, model)
    results <- rbind(results, data.table(sens_p = results_prev$sensitivity, sens_i = results_inc$sensitivity,
                                         spec_p = results_prev$specificity, spec_i = results_inc$specificity,
                                         auc_p = as.numeric(results_prev$roc.auc), auc_i = as.numeric(results_inc$roc.auc)))
  }
  results <- results[, apply(.SD, 2, mean), .SDcols = names(results)]
  table_chunk <- data.table(label = c(name, " Sensitivity", " Specificity", " AUC"),
                            `Performance for Prevalent Cases` = c("", sprintf("%.3f", results['sens_p']),
                                                             sprintf("%.3f", results['spec_p']),
                                                             sprintf("%.3f", results['auc_p'])),
                            `Performance for Incident Cases` = c("", sprintf("%.3f", results['sens_i']),
                                                                 sprintf("%.3f", results['spec_i']),
                                                                 sprintf("%.3f", results['auc_i'])))
  return(table_chunk)
}

table_dt <- rbindlist(list(get_chunk(model, "All Items"),
                      get_chunk(model2, "MMSE + I/ADLs + Word Recall"),
                      get_chunk(model3, "MMSE + I/ADLs")))

# LOOK ACROSS DIFFERENT WAVES ---------------------------------------------

get_rows <- function(model, name, prev_dt, inc_dt){
  results_prev <- get_performance(prev_dt, model)
  results_inc <- get_performance(inc_dt, model)
  result_dt <- data.table(item_set = name, sensitivity_prev = results_prev$sensitivity, sensitivity_inc = results_inc$sensitivity,
                          specificity_prev = results_prev$specificity, specificity_inc = results_inc$specificity)
  return(result_dt)
}

wave_differences <- function(wave){
  ## follow-up year
  test_dt_prev <- copy(test_dt[fu_year == wave]) ## try other follow-up years
  
  ## all incident and no-dementia cases
  test_dt_inc <- copy(test_dt[prev_dementia == 0])
  
  results <- rbind(get_rows(model, "All Items", prev_dt = test_dt_prev, inc_dt = test_dt_inc), 
                   get_rows(model2, "MMSE + I/ADLs + Word Recall", prev_dt = test_dt_prev, inc_dt = test_dt_inc),
                   get_rows(model3, "MMSE + I/ADLs", prev_dt = test_dt_prev, inc_dt = test_dt_inc))
  results[, cs_wave := wave]
  return(results)
}

allwave_results <- rbindlist(lapply(5:12, wave_differences))

# DIFFERENTIAL PERFORMANCE BY DIABETES STATUS -----------------------------

get_chunk_diabetes <- function(model, name){
  results <- data.table()
  for (wave in 5:10){
    ## follow-up year
    prev_dt <- copy(test_dt[fu_year == wave]) 
    ## all incident and no-dementia cases
    inc_dt <- copy(test_dt[prev_dementia == 0])
    results_prev_d <- get_performance(prev_dt[diabetes_tv == 1], model)
    results_inc_d <- get_performance(inc_dt[diabetes_tv == 1], model)
    results_prev_nd <- get_performance(prev_dt[diabetes_tv == 0], model)
    results_inc_nd <- get_performance(inc_dt[diabetes_tv == 0], model)
    
    results <- rbind(results, data.table(sens_p_d = results_prev_d$sensitivity, sens_i_d = results_inc_d$sensitivity,
                                         spec_p_d = results_prev_d$specificity, spec_i_d = results_inc_d$specificity,
                                         auc_p_d = as.numeric(results_prev_d$roc.auc), auc_i_d = as.numeric(results_inc_d$roc.auc),
                                         sens_p_nd = results_prev_nd$sensitivity, sens_i_nd = results_inc_nd$sensitivity,
                                         spec_p_nd = results_prev_nd$specificity, spec_i_nd = results_inc_nd$specificity,
                                         auc_p_nd = as.numeric(results_prev_nd$roc.auc), auc_i_nd = as.numeric(results_inc_nd$roc.auc)))
  }
  results <- results[, apply(.SD, 2, mean), .SDcols = names(results)]

  table_chunk <- data.table(label = c(name, " Sensitivity", " Specificity", " AUC"),
                            `Performance for Prevalent Cases - Diabetes` = c("", sprintf("%.3f", results['sens_p_d']),
                                                                             sprintf("%.3f", results['spec_p_d']),
                                                                             sprintf("%.3f", results['auc_p_d'])),
                            `Performance for Prevalent Cases - No Diabetes` = c("", sprintf("%.3f", results['sens_p_nd']),
                                                                             sprintf("%.3f", results['spec_p_nd']),
                                                                             sprintf("%.3f", results['auc_p_nd'])),
                            `Performance for Incident Cases - Diabetes` = c("", sprintf("%.3f", results['sens_i_d']),
                                                                            sprintf("%.3f", results['spec_i_d']),
                                                                            sprintf("%.3f", results['auc_i_d'])),
                            `Performance for Incident Cases - No Diabetes` = c("", sprintf("%.3f", results['sens_i_nd']),
                                                                               sprintf("%.3f", results['spec_i_nd']),
                                                                               sprintf("%.3f", results['auc_i_nd'])))
  return(table_chunk)
}

diabetes_table_dt <- rbindlist(list(get_chunk_diabetes(model, "All Items"),
                           get_chunk_diabetes(model2, "MMSE + I/ADLs + Word Recall"),
                           get_chunk_diabetes(model3, "MMSE + I/ADLs")))

# COMBINED TABLE ----------------------------------------------------------

table_dir <- paste0(dir, "paper/performance_table/")
full_table <- cbind(table_dt[,1], ## columns
                    table_dt[,2], diabetes_table_dt[,2:3], ## prevalent cases
                    table_dt[,3], diabetes_table_dt[,4:5]) ## incident cases
write.xlsx(full_table, paste0(table_dir, "full_journalpaper_table_", date, ".xlsx"))

# ROC Curves --------------------------------------------------------------

get_roc <- function(model, name, legend = F){
  
 get_prev_results <- function(wave, model_obj = model){
    print(paste0("Getting data for prev wave", wave))
    ## follow-up year
    prev_dt <- copy(test_dt[fu_year == wave]) 
    results_chunk <- get_roccuts(prev_dt, model_obj)
    results_chunk[, `:=` (measure = "prev", wave_time = wave)]
    
    return(results_chunk)
  }
  
  results_dt <- rbindlist(lapply(5:10, get_prev_results))
  results_dt[, `:=` (sens = mean(sens), spec_m1 = mean(spec_m1)), by = "cut"] ## average performance over waves
  results_dt[, wave_time := NULL]
  results_dt <- unique(results_dt)
  
  ## add incidence
  print("Getting data for inc")
  inc_dt <- copy(test_dt[prev_dementia == 0])
  results_inc <- get_roccuts(inc_dt, model)
  results_inc[, `:=` (measure = "inc")]
  results_dt <- rbind(results_dt, results_inc)
  results_dt[, label := ifelse(measure == "inc", "Performance on\nIncident Cases", "Performance on\nPrevalent Cases")]
  
  ## ROC PLOTS
  roc_plot <- ggplot(results_dt, aes(x = spec_m1, y = sens, color = label)) +
    geom_line() +
    labs(x = "1-Specificity", y = "Sensitivity") +
    scale_color_manual(values = c("Performance on\nIncident Cases" = "red",
                                  "Performance on\nPrevalent Cases" = "blue"), name = "") +
    geom_text(aes(x = 0.5, y = 0.17, label = paste0("AUC = ", auc_table[label == name, incidence])), color = "red") +
    geom_text(aes(x = 0.5, y = 0.1, label = paste0("AUC = ", auc_table[label == name, prevalence])), color = "blue") +
    ggtitle(name) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 12), legend.position = "none")
  if (legend == T){
    roc_plot <- ggplot(results_dt, aes(x = spec_m1, y = sens, color = label)) +
      geom_line() +
      labs(x = "1-Specificity", y = "Sensitivity") +
      scale_color_manual(values = c("Performance on\nIncident Cases" = "red",
                                    "Performance on\nPrevalent Cases" = "blue"), name = "") +
      geom_text(aes(x = 0.5, y = 0.17, label = paste0("AUC = ", auc_table[label == name, incidence])), color = "red") +
      geom_text(aes(x = 0.5, y = 0.1, label = paste0("AUC = ", auc_table[label == name, prevalence])), color = "blue") +
      ggtitle(name) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, size = 12), legend.position = "right")
  }  
  return(roc_plot)
}

## CREATE TABLES OF AVERAGED AUC'S TO INPUT INTO GRAPHS
auc_table <- copy(table_dt[label == " AUC" | `Performance for Prevalent Cases` == ""])
setnames(auc_table, names(auc_table), c("label", "prevalence", "incidence"))
auc_table[, `:=` (prevalence = lead(prevalence), incidence = lead(incidence))]
auc_table <- auc_table[!label == " AUC"]

roc_full <- plot_grid(get_roc(model, "All Items"), get_roc(model2, "MMSE + I/ADLs + Word Recall"),
                 get_roc(model3, "MMSE + I/ADLs", legend = T), nrow = 1, rel_widths = c(1,1,1.55))

rocplot_dir <- paste0(dir, "paper/roc_plot/")
ggsave(paste0(rocplot_dir, "main_rocs_", date, ".pdf"), plot = roc_full, width = 10, height = 3.5)


# ROC CURVES - DIABETES STATUS INCLUDED -----------------------------------

get_roc_diabetes <- function(model, name, legend = F){
  results_inc_d <- get_performance(test_dt_inc[diabetes_tv == 1], model)
  results_inc_nd <- get_performance(test_dt_inc[diabetes_tv == 0], model)
  
  ## ROC PLOTS
  roc_plot <- ggplot() +
    geom_line(aes(x = (1-results_inc_d$roc.specificities), y = results_inc_d$roc.sensitivities, color = "navyblue")) +
    geom_line(aes(x = (1-results_inc_nd$roc.specificities), y = results_inc_nd$roc.sensitivities, color = "forestgreen")) +
    labs(x = "1-Specificity", y = "Sensitivity") +
    scale_color_identity(labels = c("Diabetes", "No Diabetes"),
                         breaks = c("navyblue", "forrestgreen"),
                         name = "") +
    geom_text(aes(x = 0.5, y = 0.17, label = paste0("AUC = ", sprintf("%.3f", as.numeric(results_inc_d$roc.auc)))), color = "navyblue") +
    geom_text(aes(x = 0.5, y = 0.1, label = paste0("AUC = ", sprintf("%.3f", as.numeric(results_inc_nd$roc.auc)))), color = "forestgreen") +
    ggtitle(name) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 12))
  if (legend == T){
    roc_plot <- roc_plot +
      scale_color_identity(labels = c("Diabetes", "No Diabetes"),
                           breaks = c("navyblue", "forestgreen"),
                           name = "", guide = "legend") +
      theme(legend.spacing = unit(2, 'cm')) +
      guides(color = guide_legend(byrow = TRUE))
  }  
  return(roc_plot)
}

roc_diabetes_full <- plot_grid(get_roc_diabetes(model, "All Items"), get_roc_diabetes(model2, "MMSE + I/ADLs + Word Recall"),
                      get_roc_diabetes(model3, "MMSE + I/ADLs", legend = T), nrow = 1, rel_widths = c(1,1,1.55))

rocplot_dir <- paste0(dir, "paper/roc_plot/")
ggsave(paste0(rocplot_dir, "main_rocs_diabets_", date, ".pdf"), plot = roc_diabetes_full, width = 10, height = 3.5)


