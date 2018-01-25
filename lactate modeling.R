#create model of non-head injury with lactate and compare residuals with aishead1
library(dplyr)
library(ggplot2)
library(modelr)


coag %>% filter(hr0_lactate < 20) %>% ggplot() + 
  geom_boxplot(mapping = aes(x = as.factor(aishead1), y = hr0_lactate))

coag %>% filter(aishead1 == 0, icu_0to6h_blood_units < 30) %>% ggplot() + 
  geom_point(mapping = aes(x = hr0_lactate, y = icu_0to6h_blood_units)) + scale_y_log10()

coag %>% filter(aishead1 > 0, hr0_lactate < 20) %>% ggplot() + 
  geom_point(mapping = aes(x = hr0_lactate, y = icu_0to6h_blood_units)) + scale_y_log10()

mod1 <- lm(hr0_lactate ~ aishead1, aishead1 > 0, data = coag)
             
grid <- coag %>%
  data_grid(aishead1) %>% 
  add_predictions(mod1)

grid

coag %>% filter(aishead1 > 0, hr0_lactate < 20) %>% ggplot(aes(aishead1)) +
  geom_point(mapping = aes(y = hr0_lactate)) + 
    geom_line(data = grid, aes(y = pred))

coag <- coag %>% 
  add_residuals(mod1)

coag %>% filter(hr0_lactate <20) %>% ggplot(aes(aishead1, resid)) +
  geom_point()

coag %>% filter(hr0_lactate < 20) %>% ggplot(aes(hr0_lactate)) +
    geom_histogram()
#since its skewed, take the log of it


coag <- coag %>% filter(hr0_lactate < 20, aishead1 > 0) %>% mutate(log_lactate = log(hr0_lactate)) 

coag %>% ggplot(aes(log_lactate)) +
    geom_histogram()
#use linear regression with normal curves generally

mod2 <- lm(log_lactate ~ aishead1, data = coag)
mean_model <- lm(log_lactate ~ 1, data = coag)
summary(mean_model)
coag <- coag %>% add_residuals(mean_model)

coag %>% ggplot(aes(resid)) + geom_histogram()

summary(mod2)

coag <- coag %>% add_residuals(mod2)

coag %>% ggplot(aes(x = resid)) + geom_histogram()

coag %>% ggplot(aes(x = age, y = resid)) +
    geom_point()

coag %>% ggplot(aes(x = aishead1, y = resid)) + geom_point()
#try adding other variables - map other variables against residuals to figure out which variables to add to model

coag %>% ggplot(aes(x = aishead1, y = resid)) + geom_point() #promising, fair positive correlation
coag %>% ggplot(aes(x = hr0_resprate, y = resid)) + geom_point()
coag %>% ggplot(aes(x = hr0_sbp, y = resid)) + geom_point()
coag %>% ggplot(aes(x = hr0_factorii, y = resid)) + geom_point() #slight negative correlation
coag %>% filter(hr0_bun < 50) %>% ggplot(aes(x = hr0_bun, y = resid)) + geom_point()

mod3 <- lm(log_lactate ~ aishead1 + icu_0to6h_blood_units, data = coag)

summary(mod3) #aishead1 is a good predictor of log_lactate

coag <- coag %>% add_residuals(mod3)

coag %>% ggplot(aes(x = icu_0to6h_blood_units, y = resid)) + geom_point()

coag %>% ggplot(aes(x = hr0_factorii, y = resid)) + geom_point()
mod4 <- lm(log_lactate ~ aishead1 + icu_0to6h_blood_units + hr0_factorii, data = coag)

summary(mod4)

coag <- coag%>% add_residuals(mod4)
coag %>% ggplot(aes(x = hr0_factorii, y = resid)) + geom_point()
coag %>% ggplot(aes(x = hr0_factorii, y = resid)) + geom_point()

coag %>% filter(icu_0to6h_blood_units > 0) %>% ggplot() + 
  geom_point(mapping = aes(x = bmi, y = icu_0to6h_blood_units))
coag %>% ggplot() + 
  geom_point(mapping = aes(x = hr0_lactate, y = icu_0to6h_blood_units))
coag %>% ggplot() + geom_point(mapping = aes(x = age, y = icu_0to6h_blood_units))

bloodmod <- lm(icu_0to6h_blood_units ~ bmi, data = coag)
summary(bloodmod)

coag %>% ggplot(aes(x = age, y = icu_0to6h_blood_units)) + 
  geom_point() + stat_smooth(method = "lm", col = "red")

coag <- coag %>% add_residuals(bloodmod)
coag %>% ggplot(aes(x = age, y = resid)) + geom_point()

coag %>% filter(icu_0to6h_blood_units > 0) %>% ggplot(aes(icu_0to6h_blood_units)) + geom_histogram()
#this is a bit skewed so I'll try taking the log of it
coag <- coag %>% filter(icu_0to6h_blood_units > 0) %>% mutate(log_blood = log(icu_0to6h_blood_units))  #creates log
_blood

coag %>% ggplot(aes(log_blood)) + geom_histogram()
#this is a little bit better

log_blood_mod <- lm(log_blood ~ bmi, data = coag)
summary(log_blood_mod)

coag %>% ggplot(aes(x = bmi, y = log_blood)) + 
  geom_point() + stat_smooth(method = "lm", col = "red")

#look for papers saying certain variables are related and see if we can reproduce it using our data
#predict blood transfusion in hr0-6 - examine residuals with histograms, see what trends are apparent

blood_lactate_mod <- lm(icu_0to6h_blood_units ~ aishead1, data = coag)
summary(blood_lactate_mod)

coag %>% ggplot(aes(x=age, y = resid)) + geom_point()

#NEW

coag %>% ggplot(aes(x = hr0_lactate, y = icu_0to6h_blood_units)) + geom_point() +
    scale_y_log10()

mod5 <- lm(icu_0to6h_blood_units ~ log_lactate, data = coag)
summary(mod5) #not a very good r^2 at all

coag %>% ggplot(aes(x = bmi, y = icu_0to6h_blood_units)) + geom_point() + scale_y_log10()
coag %>% mutate(log_blood = log(icu_0to6h_blood_units)) 

mod5 <- lm(log_blood ~ hr0_lactate, data = coag)
summary(mod5) # r^2 = 0.01662

mod5 <- lm(log_blood ~ bmi, data = coag)
summary(mod5) #r^2 = 0.01498

coag %>% ggplot(aes(x = age, y = icu_0to6h_blood_units)) + geom_point() + scale_y_log10() #no real relationship

coag %>% filter(minutesToBloodDraw < 40) %>% ggplot(aes(x = minutesToBloodDraw, y = icu_0to6h_blood_units)) + 
  geom_point() + scale_y_log10()

bdmod <- lm(log_blood ~ minutesToBloodDraw, data = coag)
summary(bdmod) #r^2 = 0.153
coag %>% filter(minutesToBloodDraw < 40) %>% ggplot(aes(x = minutesToBloodDraw, y = log_blood)) + geom_point() +  
    stat_smooth(method = "lm", col = "red")

coag <- coag %>% add_residuals(bdmod)

coag %>% ggplot(aes(x = resid)) + geom_histogram() # skewed residuals

coag %>% ggplot(aes(x = bmi, y = resid)) + geom_point()

bdmod <- lm(log_blood ~ minutesToBloodDraw + hr0_lactate, data = coag)
summary(bdmod)

coag <- coag %>% add_residuals(bdmod)
coag %>% ggplot(aes(x = resid)) + geom_histogram() # still very skewed

modpen <- lm(log_blood ~ blunt, data = coag)
coag <- coag %>% add_predictions(modpen)

coag %>% ggplot(aes(x = blunt)) + geom_point(aes(y = log_blood)) +
  geom_point(aes(y = pred), color = "red", size = 4)

mod7 <- lm(log_blood ~ minutesToBloodDraw + blunt, data = coag)
summary(mod7)

coag <- coag %>% add_residuals(mod7)

coag %>% ggplot(aes(x = resid)) + geom_histogram() #need to figure out how to fix this

coag %>% ggplot(aes(x = icu_0to6h_blood_units)) + geom_histogram()
coag %>% ggplot(aes(x = log_blood)) + geom_histogram()

hr_test <- lm(log_blood ~ hr0_lactate, data = coag)
summary(hr_test)

blunt_glm <- glm(formula = icu_0to6h_blood_units ~ blunt + minutesToBloodDraw, data = coag)
#not sure if this is how to do it or how to use family parameter
summary(blunt_glm) #horrific p-value for blunt, so we should probably use something else

plot(blunt_glm)

blood_glm <- glm(formula = icu_0to6h_blood_units ~ minutesToBloodDraw + )

