##########################################################################
### Author: Emma Nichols
### Date: 1/14/2022
### Purpose: LASSO models
##########################################################################

rm(list = ls())

pacman::p_load(data.table, openxlsx, haven, readr, dplyr, magrittr, stringr,
               labelled, scales, gridExtra, grid, ggplot2, mice, survival,
               survminer, arules, gglasso)
date <- gsub("-", "_", Sys.Date())
set.seed(12494)

short <- F ## short = T will give alg for only a subset of items
sorter <- F ## shorter = T will give alg for smallest subset of items

# SET OBJECTS -------------------------------------------------------------

dir <- FILEPATH
rawdata_dir <- paste0(dir, "raw_data/")
derived_dir <- paste0(dir, "derived_data/")
model_dir <- paste0(FILEPATH, date)
plot_dir <- paste0(dir, "plots/")

# GET DATA ----------------------------------------------------------------

dt <- readr::read_rds(paste0(derived_dir, "rosmap_formatted.rds"))

variable_map <- as.data.table(read.xlsx(paste0(dir, "variable_map.xlsx")))
alg_items <- variable_map[, new_name]

if (short == T){
  alg_items <- alg_items[grepl("iadl|adl|mmse|wordrecall", alg_items)]
}
if (sorter == T){
  alg_items <- alg_items[grepl("iadl|adl|mmse", alg_items)]
}

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

# SELECT CROSS-SECTION WITH MOST CASES ------------------------------------

numcase_dt <- dt[, .(sumcase = sum(dementia, na.rm = T)), by = "fu_year"]
fu_select <- numcase_dt[sumcase == max(sumcase), fu_year]
dt <- dt[fu_year == 7 & !is.na(dementia)]

## get sample size
dt[, length(unique(projid))]

# SET UP ALGORITHM CROSS-VALIDATION ---------------------------------------

lasso_setup <- function(data){
  preds <- sort(names(data)[grepl(paste0(alg_items, collapse = "|"), names(data))])
  preds <- c(otheralg_items, preds)
  
  ## get vector of group for group lasso
  group_vector <- c(rep(1,8),2,rep(3,4),4,4,5) ## preset for demographics
  index <- 6
  for (item in sort(alg_items)){
    n <- length(preds[grepl(item, preds)])
    group_vector <- c(group_vector, rep(index, n))
    index <- index + 1
  }
  
  ## get pred matrix
  pred_mat <- as.matrix(data[, c(preds), with = F])
  
  ## get y vector
  y_vector <- data[, ifelse(dementia == 1, 1, -1)]
  
  return(list(y = y_vector, preds = pred_mat, group = group_vector))
}

# CROSS-VALIDATED GROUP LASSO ---------------------------------------------

inputs <- lasso_setup(dt)

lasso <- cv.gglasso(x = inputs$preds, y = inputs$y, group = inputs$group, pred.loss = "loss",
                        loss = "logit")

## SAVE RESULTS
dir.create(model_dir)
write_rds(lasso, paste0(model_dir, "/model.rds"))
