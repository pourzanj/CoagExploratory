library(dplyr)
library(lubridate)

#split data in to a general tbl and a tbl containing all factors and other blood tests
coag <- read.csv("Data/ACITdataset_UCSB_23Feb15.csv") %>%
  tbl_df %>%
  select(acitnum, male, age, bmi,
         mechtype, blunt, tbi, edarrivalgcs, numribfxs, admitday_intubated,
         aishead1, aisface2, aischest3, aisabdomen4, aisextremity5,
         
         prehosp_crystalloids, icu_0to6h_blood_units,
         
         injurydatetime, edarrivaldatetime, hr0_datetime, datetimeofdeath,
         
         hr0_temp, hr0_hr, hr0_sbp, hr0_basedefexc, hr0_inr,
         
         hr0_factorii, hr0_factorv, hr0_factorvii, hr0_factorviii,
         hr0_factorix, hr0_factorx, hr0_ddimer, hr0_apc, hr0_pc,
         hr0_atiii, hr0_fibrinogen_rlab) %>%
  
  #clean some of the columns
  mutate(mechtype = as.character(mechtype)) %>%
  mutate(mechtype = ifelse(mechtype == "", as.character(NA), mechtype)) %>%
  mutate(mechtype = as.factor(mechtype)) %>%
  
  mutate(blunt = as.character(blunt)) %>%
  mutate(blunt = ifelse(blunt == "", as.character(NA), blunt)) %>%
  mutate(blunt = as.factor(blunt)) %>%
  
  #wrangle times
  mutate_each(funs(parse_date_time(., "m/d/y H:M")),
              injurydatetime, edarrivaldatetime, hr0_datetime, datetimeofdeath) %>%
  
  mutate(minutesToEdArrival = (injurydatetime %--% edarrivaldatetime) / dminutes(1)) %>%
  
  #fix cases where minutes to arrival doesn't make sense because its negative
  mutate(minutesToEdArrival = ifelse(minutesToEdArrival <= 0, NA, minutesToEdArrival)) %>%
  mutate(minutesToBloodDraw = (edarrivaldatetime %--% hr0_datetime) / dminutes(1)) %>%
  mutate(minutesToBloodDraw = ifelse(!between(minutesToBloodDraw, -50, 90), NA, minutesToBloodDraw)) %>%
  mutate(daysToDeath = (edarrivaldatetime %--% datetimeofdeath) / ddays(1)) %>%
  mutate(died = !is.na(datetimeofdeath)) %>%
  select(-injurydatetime, -hr0_datetime, -edarrivaldatetime, -datetimeofdeath)