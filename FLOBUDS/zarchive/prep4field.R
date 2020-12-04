Intercept <- coef(modelhere, prob=c(0.25, 0.75))$species[, c(1, 3:4), 1] %>% ## here we make find the posterior distributions and means for each predictor
  as.data.frame() %>%
  round(digits = 2) %>% 
  rename(mean = Estimate) %>%
  rename(`25%` = Q25) %>%
  rename(`75%` = Q75) %>%
  dplyr::select(mean, `25%`, `75%`) ### can change according to uncertainty intervals you want
new.names<-NULL
for(i in 1:length(spp)){
  new.names[i]<-paste("Intercept")}

year1990$parameter<-new.names
year1991 <- coef(modelhere, prob=c(0.25, 0.75))$species[, c(1, 3:4), 2] %>%
  as.data.frame() %>%
  round(digits = 2) %>% 
  rename(mean = Estimate) %>%
  rename(`25%` = Q25) %>%
  rename(`75%` = Q75) %>%
  dplyr::select( mean, `25%`, `75%`) 
new.names<-NULL
for(i in 1:length(spp)){
  new.names[i]<-paste("1991")
}
year1991$parameter<-new.names
mod.ranef<-full_join(year1990,year1991)

year1992 <- coef(modelhere, prob=c(0.25, 0.75))$species[, c(1, 3:4), 3] %>%
  as.data.frame() %>%
  round(digits = 2) %>% 
  rename(mean = Estimate) %>%
  rename(`25%` = Q25) %>%
  rename(`75%` = Q75) %>%
  dplyr::select( mean, `25%`, `75%`) 
new.names<-NULL
for(i in 1:length(spp)){
  new.names[i]<-paste("1992")
}

year1992$parameter<-new.names
mod.ranef <- full_join(mod.ranef, year1992)

year1993 <- coef(modelhere, prob=c(0.25, 0.75))$species[, c(1, 3:4), 4] %>%
  as.data.frame() %>%
  round(digits = 2) %>% 
  rename(mean = Estimate) %>%
  rename(`25%` = Q25) %>%
  rename(`75%` = Q75) %>%
  dplyr::select( mean, `25%`, `75%`) 
new.names<-NULL
for(i in 1:length(spp)){
  new.names[i]<-paste("1993")
}

year1993$parameter<-new.names
mod.ranef<-full_join(mod.ranef, year1993)

year1994 <- coef(modelhere, prob=c(0.25, 0.75))$species[, c(1, 3:4), 5] %>%
  as.data.frame() %>%
  round(digits = 2) %>% 
  rename(mean = Estimate) %>%
  rename(`25%` = Q25) %>%
  rename(`75%` = Q75) %>%
  dplyr::select( mean, `25%`, `75%`) 
new.names<-NULL
for(i in 1:length(spp)){
  new.names[i]<-paste("1994")
}



year1994$parameter<-new.names
mod.ranef <- full_join(mod.ranef, year1994)

year1995 <- coef(modelhere, prob=c(0.25, 0.75))$species[, c(1, 3:4), 6] %>%
  as.data.frame() %>%
  round(digits = 2) %>% 
  rename(mean = Estimate) %>%
  rename(`25%` = Q25) %>%
  rename(`75%` = Q75) %>%
  dplyr::select( mean, `25%`, `75%`) 
new.names<-NULL
for(i in 1:length(spp)){
  new.names[i]<-paste("1995")
}

year1995$parameter<-new.names
mod.ranef <- full_join(mod.ranef, year1995)

year1996 <- coef(modelhere, prob=c(0.25, 0.75))$species[, c(1, 3:4), 7] %>%
  as.data.frame() %>%
  round(digits = 2) %>% 
  rename(mean = Estimate) %>%
  rename(`25%` = Q25) %>%
  rename(`75%` = Q75) %>%
  dplyr::select( mean, `25%`, `75%`) 
new.names<-NULL
for(i in 1:length(spp)){
  new.names[i]<-paste("1996")
}

year1996$parameter<-new.names
mod.ranef <- full_join(mod.ranef, year1996)


###############################################



library(broom)
modoutput1 <- tidy(modelhere, prob=c(0.5))

