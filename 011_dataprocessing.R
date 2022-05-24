##########################################################################
### Author: Emma Nichols
### Date: 8/12/2021
### Purpose: ROSMAP Data Prep
##########################################################################

rm(list = ls())

pacman::p_load(data.table, openxlsx, haven, readr, dplyr, magrittr, stringr,
               labelled, scales, gridExtra, grid, ggplot2)
date <- gsub("-", "_", Sys.Date())
set.seed(6541)

# SET OBJECTS -------------------------------------------------------------

dir <- FILEPATH
rawdata_dir <- paste0(dir, "raw_data/")
derived_dir <- paste0(dir, "derived_data/")

# GET DATA ----------------------------------------------------------------

baseline_dt <- as.data.table(read.xlsx(paste0(rawdata_dir, "dataset_1014_basic_08-20-2021.xlsx")))
long_dt <- as.data.table(read.xlsx(paste0(rawdata_dir, "dataset_1014_long_08-20-2021.xlsx")))

rosmap_dt <- merge(baseline_dt, long_dt, all = T, by = c("projid", "study"))

variable_map <- as.data.table(read.xlsx(paste0(dir, "variable_map.xlsx")))

# MISSING DATA ------------------------------------------------------------

## FILL IN MISSING DIAGNOSIS IF HAD SAME DIAGNOSIS AT VISIT BEFORE AND AFTER
rosmap_dt[, rowid := 1:.N]
ids <- rosmap_dt[is.na(dcfdx), unique(rowid)]
for (id in ids){
  pid <- rosmap_dt[rowid == id, projid]
  fu <- rosmap_dt[rowid == id, fu_year]
  before_dx <- rosmap_dt[fu_year == (fu-1) & projid == pid, dcfdx]
  after_dx <- rosmap_dt[fu_year == (fu+1) & projid == pid, dcfdx]
  if (length(before_dx) == 1 & length(after_dx) == 1){
    if (!is.na(before_dx) & !is.na(after_dx)){
      if (before_dx == after_dx) {
        rosmap_dt[rowid == id, dcfdx := before_dx]
      }   
    }
  }
}

# GET ITEM-LEVEL DATA -----------------------------------------------------

# COGNIGITVE DATA ---------------------------------------------------------

ros_cog <- as.data.table(read_delim(paste0(rawdata_dir, "ros/cts_update.txt"), 
                                    delim = "|"))
map_cog <- as.data.table(read_delim(paste0(rawdata_dir, "map/cts_update.txt"),
                                    delim = "|"))

cog_dt <- rbind(ros_cog, map_cog)

## FIX MISSINGNESS FOR DIGIT SPAN FORWARDS
## TESTING STOPS AFTER TWO CONSECUTIVE ERRORS AT GIVEN SEQUENCE LENGTH
fix_digit <- function(dt, items){
  dt[, (items) := lapply(.SD, function(x) ifelse(x == 9, 0, x)), .SDcols = items]
  for (item_round in 1:5){
    round_items <- items[grepl(paste0(item_round, "$"), items)]
    dt[, paste0("round", item_round) := rowSums(.SD), .SDcols = round_items]
    change_items <- items[grepl(paste0("[", item_round + 1, "-6]$"), items)]
    dt[get(paste0("round", item_round)) == 0, (change_items) := lapply(.SD, function(x) ifelse(is.na(x), 0, x)), .SDcols = change_items]
  }
  dt[, paste0("round", 1:5) := NULL]
  return(dt)
}

cog_dt <- fix_digit(cog_dt, c(paste0("digfora_item", 1:6), paste0("digforb_item", 1:6)))
cog_dt <- fix_digit(cog_dt, c(paste0("digbaka_item", 1:6), paste0("digbakb_item", 1:6)))

## FIX MISSINGNESS FOR NUMBER COMPARISON
## IF MISSING ON WRONG VARIBALES - THIS IS TURNED INTO "WRONG"

c_witems <- paste0("nc_wrong_item", 1:6)
cog_dt[, (c_witems) := lapply(.SD, function(x) ifelse(is.na(x), 1, x)), .SDcols = c_witems]

## FIX MISSINGNESS FOR DIGIT ORDER
dt <- copy(cog_dt)

fix_digitorder <- function(dt){
  do_items <- c(paste0("do_item", 1:14))
  dt[, summissing := apply(.SD, 1, function(x) sum(x == 8899, na.rm = T)), .SDcols = do_items]
  dt[summissing == 0, (do_items) := lapply(.SD, function(x) ifelse(is.na(x), 0, x)), .SDcols = do_items]
  dt[, (do_items) := lapply(.SD, function(x) ifelse(x == 8899, NA, x)), .SDcols = do_items]
  dt[, c("summissing") := NULL]
}

cog_dt <- fix_digitorder(cog_dt)

## FIX CATEGORY FLUENCY
cog_dt[is.na(animals_item1) & !is.na(fruits_item1), animals_item1 := fruits_item1]
cog_dt[!is.na(animals_item1) & is.na(fruits_item1), fruits_item1 := animals_item1]

# PROCESS AND RENAME ALL DATA ---------------------------------------------

process_data <- function(dt, f){
  map <- variable_map[file == f]
  add_map <- map[grepl("\\+|\\-", orig_name)]
  if (nrow(add_map) > 0){
    for (row in 1:nrow(add_map)){
      new <- add_map[, new_name][row]
      equation <- add_map[, orig_name][row]
      dt[, c(new) := eval(parse(text = equation))]
    }
  }
  other_map <- map[!grepl("\\+|\\-", orig_name)]
  for (row in 1:nrow(other_map)){
    new <- other_map[, new_name][row]
    old <- other_map[, orig_name][row]
    dt[, c(new) := get(old)]
  }
  dt <- dt[, c("projid", "fu_year", map[, new_name]), with = F]
  return(dt)
}

## COGNITION
cog_dt <- process_data(cog_dt, f = "cts")

## ADLS
ros_adl <- as.data.table(read_delim(paste0(rawdata_dir, "ros/katz.txt"), 
                                    delim = "|"))
map_adl <- as.data.table(read_delim(paste0(rawdata_dir, "map/katz.txt"),
                                    delim = "|"))
adl_dt <- rbind(ros_adl, map_adl)
adl_dt <- process_data(adl_dt, f = "katz")

##IADLS
ros_iadl <- as.data.table(read_delim(paste0(rawdata_dir, "ros/iadl.txt"), 
                                    delim = "|"))
map_iadl <- as.data.table(read_delim(paste0(rawdata_dir, "map/iadl.txt"),
                                    delim = "|"))
iadl_dt <- rbind(ros_iadl, map_iadl)
iadl_dt <- process_data(iadl_dt, f = "iadl")

# DEMOGRAPHIC VARIABLES ---------------------------------------------------

## DEMENTIA STATUS
rosmap_dt[, dementia := as.numeric(dcfdx %in% c(4:6))]
rosmap_dt[is.na(dcfdx), dementia := NA]

## COLLAPSE RACE
race_dt <- data.table(race7 = c(1:7, NA),
                      race = factor(c("white", "black", rep("other", 5), NA),
                                    levels = c("white", "black", "other")))
rosmap_dt <- merge(rosmap_dt, race_dt, by = "race7", all.x = T)

## HISPANIC
rosmap_dt[, hispanic := as.numeric(spanish == 1)]

# DIABETES AND RISKS ------------------------------------------------------

## smoking 
rosmap_dt[, smoke_cf := as.numeric(!smoking == 0)]

## baseline hypertension 
rosmap_dt[, hyp_bl := hypertension_bl]

## claudication 
rosmap_dt[, claud_bl := claudication_bl]

## heart attack baseline 
rosmap_dt[, heart_bl := heart_bl]

## stroke baseline
rosmap_dt[, stroke_bl := stroke_bl]

## rename diabetes - full calculation below
rosmap_dt[, diabetes_tv := dm_cum]

## baseline bmi 
bmi_dt <- copy(rosmap_dt[fu_year == 0, .(projid, bmi)])
setnames(bmi_dt, "bmi", "bmi_bl")
rosmap_dt <- merge(rosmap_dt, bmi_dt, by = "projid", all.x = T)

rosmap_dt <- rosmap_dt[, .(projid, fu_year, study, age_at_visit, age_bl, age_death, educ, msex, race, hispanic, dementia, cogn_global,
                           smoke_cf, hyp_bl, claud_bl, heart_bl, stroke_bl, bmi_bl, diabetes_tv, hba1c, mmse = cts_mmse30)]

# PUT IT ALL TOGETHER -----------------------------------------------------

rosmap_dt <- Reduce(function(x,y) merge(x, y, all.x = T, by = c("projid", "fu_year")), list(rosmap_dt, cog_dt, adl_dt, iadl_dt))

# DEFINE DIABETES ---------------------------------------------------------

## define diabetes (combine hba1c with diabetes_tv)
rosmap_dt[, hba1c_pos := as.numeric(hba1c > 6.5)]
rosmap_dt[diabetes_tv == 0 & !is.na(hba1c_pos), diabetes_tv := ifelse(hba1c_pos == 1, 1, diabetes_tv)] 
suppressWarnings(rosmap_dt[, any_diabetes := as.numeric(max(diabetes_tv, na.rm = T) == 1), by = "projid"])
rosmap_dt[diabetes_tv == 1, first_diabetes := min(fu_year), by = "projid"]
rosmap_dt[any_diabetes == 1, first_diabetes := max(first_diabetes, na.rm = T), by = "projid"] ## create variable for all visits
rosmap_dt[is.na(first_diabetes), first_diabetes := 100] ## impossibly large indicator for those without diabetes
## fix people who revert (dim(dt[fu_year > first_diabetes & diabetes_tv == 0]))
rosmap_dt[any_diabetes == 1 & fu_year >= first_diabetes, diabetes_tv := 1]

# DEFINE INCIDENT AND PREVALENT DEMENTIA ----------------------------------

rosmap_dt <- rosmap_dt[order(projid, fu_year)]

## INCIDENT DEMENTIA
suppressWarnings(rosmap_dt[, any_dementia := as.numeric(max(dementia, na.rm = T)), by = "projid"])
rosmap_dt[is.infinite(any_dementia), any_dementia := 0]
rosmap_dt[dementia == 1, first_dem_visit := min(fu_year), by = "projid"]
rosmap_dt[any_dementia == 1, first_dem_visit := max(first_dem_visit, na.rm = T), by = "projid"] ## set first dementia visit for all visits among those with dementia

## GET RID OF INDIVIDUALS WITH PREVALENT DEMENTIA AT FIRST VISIT
message(paste0("getting rid of ", nrow(rosmap_dt[first_dem_visit == 0]), " rows, ",
               rosmap_dt[first_dem_visit == 0, length(unique(projid))]," people, with prevalent dementia at baseline"))
rosmap_dt <- rosmap_dt[first_dem_visit > 0 | is.na(first_dem_visit)]

rosmap_dt[is.na(first_dem_visit), first_dem_visit := 100] ## set to arbitrarily large and non-existent visit
rosmap_dt[, inc_dementia := as.numeric(first_dem_visit == fu_year)]

## COUNT NUMBER OF TIMES WE SEE REVERTING FROM DEMENTIA 
message(paste0("Number of rows that revert from dementia: ", nrow(rosmap_dt[any_dementia == 1 & fu_year >= first_dem_visit & !dementia == 1]),
               " (", sprintf("%.2f", nrow(rosmap_dt[any_dementia == 1 & fu_year >= first_dem_visit & !dementia == 1])/nrow(rosmap_dt)*100), "% of total rows)", 
               ", Number of individuals that revert from dementia: ", rosmap_dt[any_dementia == 1 & fu_year >= first_dem_visit & !dementia == 1, length(unique(projid))],
               " (", sprintf("%.2f", rosmap_dt[any_dementia == 1 & fu_year >= first_dem_visit & !dementia == 1, length(unique(projid))]/rosmap_dt[, length(unique(projid))]*100),
               "% of individuals in the study)"))

## MAKE SURE ALL VISITS ARE DEMENTIA AFTER FIRST DEMENTIA VISIT
rosmap_dt[any_dementia == 1 & fu_year >= first_dem_visit, dementia := 1]

## PREVALENT DEMENTIA
rosmap_dt[, prev_dementia := as.numeric(dementia == 1 & ! inc_dementia == 1)]

# EXCLUDE PEOPLE WITH MISSING DEMOGRAPHICS --------------------------------

nrows_missingage <- nrow(rosmap_dt[is.na(age_at_visit)])
nid_missingage <- rosmap_dt[is.na(age_at_visit), length(unique(projid))]
message(paste0("Dropping ", nrows_missingage, " visits from ", nid_missingage, " people for missing age at visit"))
rosmap_dt <- rosmap_dt[!is.na(age_at_visit)]

nrows_missingeduc <- nrow(rosmap_dt[is.na(educ)])
nid_missingeduc <- rosmap_dt[is.na(educ), length(unique(projid))]
message(paste0("Dropping ", nrows_missingeduc, " visits from ", nid_missingeduc, " people for missing education at visit"))
rosmap_dt <- rosmap_dt[!is.na(educ)]

# DEFINE TRAIN AND TEST SETS ----------------------------------------------

## try selecting 1/3 of individuals
projids <- rosmap_dt[, unique(projid)]
selected_projids <- sample(x = projids, size = round(length(projids)/3), replace = F)
rosmap_dt <- rosmap_dt[, testdata := ifelse(projid %in% selected_projids, 1, 0)]

# SAVE INTERIM DATA -------------------------------------------------------

readr::write_rds(rosmap_dt, paste0(derived_dir, "rosmap_formatted.rds"))
write.csv(rosmap_dt, paste0(derived_dir, "rosmap_formatted.csv"), row.names = F)
