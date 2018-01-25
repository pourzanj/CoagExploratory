ggplot(coag) + geom_histogram(aes(x = icu_0to6h_blood_units))

ggplot(coag) + geom_boxplot(aes(x = blunt, y= icu_0to6h_blood_units))

coag %>% 
  group_by(bmi) %>% 
  summarise(m = mean(icu_0to6h_blood_units, na.rm = TRUE), 
            ct = n(), 
            na = sum(is.na(icu_0to6h_blood_units)))
