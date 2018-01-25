library(dplyr)
library(ggplot2)
library(modelr)


pen <- filter(coag, blunt == "Penetrating")
print(pen)


BmiHist <- ggplot(data = coag) + geom_histogram(mapping = aes(x = bmi))
BmiHist + scale_x_log10()
#coag$bmigr<-cut(coag$bmi, c(0,5,10,15,20,25,30,35,40,45,50))

coag <- coag %>% mutate(bmigr = cut(bmi, c(0,5,10,15,20,25,30,35,40,45,50)))

coag <- coag %>% mutate(bmigr2 = cut(bmi, quantile(bmi, na.rm = TRUE)))


col <- select(pen, age, bmi, bmigr, mechtype, icu_0to6h_blood_units, blunt, numribfxs, daysToDeath, died)

BloodVsBmi <- ggplot(data = col) + geom_point(mapping = aes(x = bmi, y = icu_0to6h_blood_units))
BloodVsBmi + ylim(0,20)
BloodVsBmi + scale_y_log10()


ggplot(data = col) + geom_bin2d(mapping = aes(x = bmi, y = icu_0to6h_blood_units))


ggplot(data = coag) + geom_bar(mapping = aes(x=bmigr))

quantile(coag$bmi, na.rm = TRUE)

coag %>% 
  group_by(bmigr) %>%
  summarise(m = mean(icu_0to6h_blood_units, na.rm = TRUE), ct = n(), na = sum(is.na(icu_0to6h_blood_units)))
ggplot(data = coag) + geom_point(mapping = aes(x = bmigr, y = mean(icu_0to6h_blood_units))) #why doesn't this work?

ggplot(data = coag) + geom_boxplot(aes(x = bmigr, y = icu_0to6h_blood_units))

ggplot(data = coag) + geom_boxplot(aes (x = bmigr2, y = icu_0to6h_blood_units)) +
  scale_y_log10()

coag %>% 
    group_by(mechtype) %>%
        summarise(m = mean(icu_0to6h_blood_units, na.rm = TRUE), ct = n(), na = sum(is.na(icu_0to6h_blood_units)))


ggplot(data = coag) + geom_bar(aes(x = mechtype))
ggplot(data = coag) + geom_boxplot(aes(x = mechtype, y = icu_0to6h_blood_units)) +
  scale_y_log10() +
  facet_grid(died ~ tbi)

# Different Ways to Visualize Multiple Variables
# 1) Colors and Shapes
# 2) Faceting
# 3) Group By (possibly multiple variables) and look at means std dev.

coag %>% 
  group_by(died, bmigr2) %>%
  summarise(mean_dead_blood = mean(icu_0to6h_blood_units, na.rm = TRUE),
            mean_tbi = mean(tbi),
            mean_age = mean(age),
            mean_head_injury = mean(aishead1),
            mean_bp = mean(hr0_sbp),
            mean_basedef = mean(hr0_basedefexc),
            ct = n(),
            na = sum(is.na(icu_0to6h_blood_units)))

ggplot(data = col) + geom_point(mapping = aes(x = bmi, y = daysToDeath)) + 
    scale_y_log10()

ggplot(data = col) + geom_bar(aes(x = died))
ggplot(data = col) + geom_boxplot(aes(x = died, y = icu_0to6h_blood_units))

ggplot(data = col) +
  geom_count(mapping = aes(x = died, y = icu_0to6h_blood_units))

#grouping ages by quantiles
coag <- coag %>% mutate(agegr2 = cut(age, quantile(age, na.rm = TRUE)))

#observing instances of ages in data by death type
ggplot(data = coag) +
  geom_histogram(mapping = aes(x = age)) +
  facet_wrap(~ mechtype, nrow = 2)

#how do bmi vs. blood vary separated by death type
ggplot(data = coag) +
  geom_count(mapping = aes(x = bmi, y  = icu_0to6h_blood_units)) +
    facet_wrap(~ mechtype, nrow = 2)

#finding mean statistics based on "blunt" and age quantiles
coag %>% 
  group_by(agegr2) %>%
  summarise(mean_blood = mean(icu_0to6h_blood_units, na.rm = TRUE),
            mean_head_injury = mean(aishead1),
            mean_bp = mean(hr0_sbp),
            mean_basedef = mean(hr0_basedefexc),
            ct = n(),
            na = sum(is.na(icu_0to6h_blood_units)))


coag %>%
  filter(mechtype == "Fall/Crush") %>% group_by (agegr2) %>%
    summarise(mean_blood = mean(icu_0to6h_blood_units, na.rm = TRUE),
              mean_head_injury = mean(aishead1),
              mean_bp = mean(hr0_sbp),
              mean_basedef = mean(hr0_basedefexc),
              ct = n(),
              na = sum(is.na(icu_0to6h_blood_units)))

coag %>%
    filter(mechtype == "GSW") %>% group_by(agegr2) %>% 
      summarise(mean_blood = mean(icu_0to6h_blood_units, na.rm = TRUE),
                mean_head_injury = mean(aishead1),
                mean_bp = mean(hr0_sbp),
                mean_basedef = mean(hr0_basedefexc),
                ct = n(),
                na = sum(is.na(icu_0to6h_blood_units)))
#In fall/crush wounds - blood transferred seems to go down with age, while it increases with age in
#Gunshot wounds - also gunshot wounds have higher average units of blood transferred than fall/crush
#probably because its a penetrating wound rather than blunt
coag %>% filter(mechtype == "GSW" | mechtype == "Fall/Crush") %>%
ggplot() +
    geom_boxplot(mapping = aes(x= agegr2, y = icu_0to6h_blood_units)) +
      #scale_y_log10() + 
  facet_wrap(~ mechtype )

coag %>% filter(mechtype == "GSW" | mechtype == "Fall/Crush") %>%
    ggplot() +
    geom_boxplot(mapping = aes(x = agegr2, y = icu_0to6h_blood_units)) +
    scale_y_log10() + facet_grid(as.factor(aishead1)~mechtype)

#extremities seem to have significantly higher blood loss

coag %>% ggplot() + 
    geom_boxplot(mapping = aes(x = agegr2, y = icu_0to6h_blood_units)) + 
      + scale_y_log10()

#trying to see if any observable relationship between lactate concentration and blood given -
#read online that it should be correlated with tbi and there is slightly observable positive correlation
coag %>% ggplot() + 
    geom_point(mapping = aes(x = hr0_lactate, y = icu_0to6h_blood_units)) +
      scale_y_log10()
 
#look at how lactate varies with BP
coag %>% ggplot() + 
  geom_point(mapping = aes(x = hr0_lactate, y = icu_0to6h_blood_units)) + scale_y_log10() + 
  facet_grid(tbi~died)

coag  %>% filter(!is.na(tbi)) %>% ggplot() + 
    geom_density(mapping = aes(x = hr0_lactate, fill = as.factor(tbi)), alpha = 0.5)

#read also that gunshot wounds specifically in abdomen use lactate as measure of severity of injury.
#sample size too small to conclude anything
coag %>% filter(mechtype == "GSW")  %>%
  ggplot() + geom_boxplot(mapping = aes(x = as.factor(aisabdomen4), y = hr0_lactate))

#lactate seems to be more of a measure in GSW's than in fall/crush
coag %>% filter(mechtype == "GSW" | mechtype == "Fall/Crush") %>%
  ggplot() + geom_point(mapping = aes(x = hr0_lactate, y = icu_0to6h_blood_units)) + 
    facet_wrap(~mechtype) + scale_y_log10()

coag %>% filter(blunt == "Blunt") %>%
  ggplot() + geom_point(mapping = aes(x = hr0_lactate, y = icu_0to6h_blood_units)) + 
  scale_y_log10()

#delving deeper into lactate specifically - not necessarily related to blood yet but i will
#tie it back later: lactate seems to increase with severity of head injury
coag %>% filter(hr0_lactate < 20) %>% ggplot() + 
  geom_boxplot(mapping = aes(x = as.factor(aishead1), y = hr0_lactate))
#create model of non-head injury with lactate and compare residuals with aishead1
 #look at medians in summarise 
coag %>% group_by(as.factor(aishead1)) %>%
    summarise(mean_blood = mean(icu_0to6h_blood_units, na.rm = TRUE),
              mean_lactate = mean(hr0_lactate, na.rm = TRUE),
              mean_bp = mean(hr0_sbp, na.rm = TRUE),
              mean_basedef = mean(hr0_basedefexc, na.rm = TRUE),
              ct = n(),
              na = sum(is.na(icu_0to6h_blood_units)))
  
#lactate seems to increase with severity of wound
  
coag %>% filter(mechtype == "GSW")  %>%
  ggplot() + geom_point(mapping = aes(x = iss, y = hr0_lactate))

coag <- coag %>% mutate(lactgr = cut(hr0_lactate, c(0,2,4,6,8,10,12,14)))

coag %>% group_by(lactgr) %>%
  summarise(mean_blood = mean(icu_0to6h_blood_units, na.rm = TRUE),
            mean_bp = mean(hr0_sbp, na.rm = TRUE),
            mean_basedef = mean(hr0_basedefexc, na.rm = TRUE),
            ct = n(),
            na = sum(is.na(icu_0to6h_blood_units)))

coag %>% ggplot() + 
    geom_boxplot(mapping = aes(x = lactgr, y = icu_0to6h_blood_units)) + scale_y_log10()

coag %>% ggplot() + 
  geom_boxplot(mapping = aes(x = lactgr, y = icu_0to6h_blood_units)) + scale_y_log10() + 
  facet_wrap(~died)

#similar analysis done with age - seems pretty inconclusive
ggplot(data = coag) +
  geom_count(mapping = aes(x = age, y  = icu_0to6h_blood_units)) +
  facet_wrap(~ blunt, nrow = 2)

#I'm curious about the line of red's at 50
ggplot(data = coag) + 
    geom_count(mapping = aes(x = prehosp_crystalloids, y = icu_0to6h_blood_units, color = died)) + 
      scale_y_log10()

#Read that death may vary with crystalloids
ggplot(data = coag) + 
  geom_count(mapping = aes(x = prehosp_crystalloids, y = daysToDeath))

ggplot(data = coag) + 
    geom_histogram(mapping = aes(x = prehosp_crystalloids))

#mean statistics on blood and crystalloids by death
coag %>% 
    group_by(died) %>%
summarise(mean_dead_blood = mean(icu_0to6h_blood_units, na.rm = TRUE),
          ct = n(),
          na = sum(is.na(icu_0to6h_blood_units)))

coag %>% 
  group_by(died) %>%
  summarise(mean_dead_cryst = mean(prehosp_crystalloids, na.rm = TRUE),
            ct = n(),
            na = sum(is.na(icu_0to6h_blood_units)))


