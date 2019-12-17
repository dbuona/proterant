Intercept <- coef(modelhere, prob=c(0.25, 0.75))$GEN.SPA[, c(1, 3:4), 1] %>% ## here we make find the posterior distributions and means for each predictor
  as.data.frame() %>%
  round(digits = 2) %>% 
  rename(mean = Estimate) %>%
  rename(`25%` = Q25) %>%
  rename(`75%` = Q75) %>%
  dplyr::select(mean, `25%`, `75%`) ### can change according to uncertainty intervals you want
new.names<-NULL
for(i in 1:length(spp)){
  new.names[i]<-paste("Intercept", "[", i, "]", sep="")}

Intercept$parameter<-new.names
Chill <- coef(modelhere, prob=c(0.25, 0.75))$GEN.SPA[, c(1, 3:4), 2] %>%
  as.data.frame() %>%
  round(digits = 2) %>% 
  rename(mean = Estimate) %>%
  rename(`25%` = Q25) %>%
  rename(`75%` = Q75) %>%
  dplyr::select( mean, `25%`, `75%`) 
new.names<-NULL
for(i in 1:length(spp)){
  new.names[i]<-paste("Chill", "[", i, "]", sep="")
}
Chill$parameter<-new.names
mod.ranef<-full_join(Intercept, Chill)
Light <- coef(modelhere, prob=c(0.25, 0.75))$GEN.SPA[, c(1, 3:4), 3] %>%
  as.data.frame() %>%
  round(digits = 2) %>% 
  rename(mean = Estimate) %>%
  rename(`25%` = Q25) %>%
  rename(`75%` = Q75) %>%
  dplyr::select( mean, `25%`, `75%`) 
new.names<-NULL
for(i in 1:length(spp)){
  new.names[i]<-paste("Light", "[", i, "]", sep="")
}

Light$parameter<-new.names
mod.ranef <- full_join(mod.ranef, Light)

Force <- coef(modelhere, prob=c(0.25, 0.75))$GEN.SPA[, c(1, 3:4), 4] %>%
  as.data.frame() %>%
  round(digits = 2) %>% 
  rename(mean = Estimate) %>%
  rename(`25%` = Q25) %>%
  rename(`75%` = Q75) %>%
  dplyr::select( mean, `25%`, `75%`) 
new.names<-NULL
for(i in 1:length(spp)){
  new.names[i]<-paste("Force", "[", i, "]", sep="")
}

Force$parameter<-new.names
mod.ranef<-full_join(mod.ranef, Force)

Chill_Light <- coef(modelhere, prob=c(0.25, 0.75))$GEN.SPA[, c(1, 3:4), 5] %>%
  as.data.frame() %>%
  round(digits = 2) %>% 
  rename(mean = Estimate) %>%
  rename(`25%` = Q25) %>%
  rename(`75%` = Q75) %>%
  dplyr::select( mean, `25%`, `75%`) 
new.names<-NULL
for(i in 1:length(spp)){
  new.names[i]<-paste("Chill:Light", "[", i, "]", sep="")
}



Chill_Light$parameter<-new.names
mod.ranef <- full_join(mod.ranef, Chill_Light)

Chill_Force <- coef(modelhere, prob=c(0.25, 0.75))$GEN.SPA[, c(1, 3:4), 6] %>%
  as.data.frame() %>%
  round(digits = 2) %>% 
  rename(mean = Estimate) %>%
  rename(`25%` = Q25) %>%
  rename(`75%` = Q75) %>%
  dplyr::select( mean, `25%`, `75%`) 
new.names<-NULL
for(i in 1:length(spp)){
  new.names[i]<-paste("Chill:Force", "[", i, "]", sep="")
}

Chill_Force$parameter<-new.names
mod.ranef <- full_join(mod.ranef, Chill_Force)

Light_Force <- coef(modelhere, prob=c(0.25, 0.75))$GEN.SPA[, c(1, 3:4), 7] %>%
  as.data.frame() %>%
  round(digits = 2) %>% 
  rename(mean = Estimate) %>%
  rename(`25%` = Q25) %>%
  rename(`75%` = Q75) %>%
  dplyr::select( mean, `25%`, `75%`) 
new.names<-NULL
for(i in 1:length(spp)){
  new.names[i]<-paste("Light:Force", "[", i, "]", sep="")
}

Light_Force$parameter<-new.names
mod.ranef <- full_join(mod.ranef, Light_Force)


###############################################


Intercept <- coef(modelhere2, prob=c(0.25, 0.75))$GEN.SPA[, c(1, 3:4), 1] %>% ## here we make find the posterior distributions and means for each predictor
  as.data.frame() %>%
  round(digits = 2) %>% 
  rename(mean = Estimate) %>%
  rename(`25%` = Q25) %>%
  rename(`75%` = Q75) %>%
  dplyr::select(mean, `25%`, `75%`) ### can change according to uncertainty intervals you want
new.names<-NULL
for(i in 1:length(spp)){
  new.names[i]<-paste("Intercept", "[", i, "]", sep="")}

Intercept$parameter<-new.names
Chill <- coef(modelhere2, prob=c(0.25, 0.75))$GEN.SPA[, c(1, 3:4), 2] %>%
  as.data.frame() %>%
  round(digits = 2) %>% 
  rename(mean = Estimate) %>%
  rename(`25%` = Q25) %>%
  rename(`75%` = Q75) %>%
  dplyr::select( mean, `25%`, `75%`) 
new.names<-NULL
for(i in 1:length(spp)){
  new.names[i]<-paste("Chill", "[", i, "]", sep="")
}
Chill$parameter<-new.names
mod.ranef2<-full_join(Intercept, Chill)
Light <- coef(modelhere2, prob=c(0.25, 0.75))$GEN.SPA[, c(1, 3:4), 3] %>%
  as.data.frame() %>%
  round(digits = 2) %>% 
  rename(mean = Estimate) %>%
  rename(`25%` = Q25) %>%
  rename(`75%` = Q75) %>%
  dplyr::select( mean, `25%`, `75%`) 
new.names<-NULL
for(i in 1:length(spp)){
  new.names[i]<-paste("Light", "[", i, "]", sep="")
}

Light$parameter<-new.names
mod.ranef2 <- full_join(mod.ranef2, Light)
Force <- coef(modelhere2, prob=c(0.25, 0.75))$GEN.SPA[, c(1, 3:4), 4] %>%
  as.data.frame() %>%
  round(digits = 2) %>% 
  rename(mean = Estimate) %>%
  rename(`25%` = Q25) %>%
  rename(`75%` = Q75) %>%
  dplyr::select( mean, `25%`, `75%`) 
new.names<-NULL
for(i in 1:length(spp)){
  new.names[i]<-paste("Force", "[", i, "]", sep="")
}

Force$parameter<-new.names
mod.ranef2<-full_join(mod.ranef2, Force)

Chill_Light <- coef(modelhere2, prob=c(0.25, 0.75))$GEN.SPA[, c(1, 3:4), 5] %>%
  as.data.frame() %>%
  round(digits = 2) %>% 
  rename(mean = Estimate) %>%
  rename(`25%` = Q25) %>%
  rename(`75%` = Q75) %>%
  dplyr::select( mean, `25%`, `75%`) 
new.names<-NULL
for(i in 1:length(spp)){
  new.names[i]<-paste("Chill:Light", "[", i, "]", sep="")
}

Chill_Light$parameter<-new.names
mod.ranef2 <- full_join(mod.ranef2, Chill_Light)

Chill_Force <- coef(modelhere2, prob=c(0.25, 0.75))$GEN.SPA[, c(1, 3:4), 6] %>%
  as.data.frame() %>%
  round(digits = 2) %>% 
  rename(mean = Estimate) %>%
  rename(`25%` = Q25) %>%
  rename(`75%` = Q75) %>%
  dplyr::select( mean, `25%`, `75%`) 
new.names<-NULL
for(i in 1:length(spp)){
  new.names[i]<-paste("Chill:Force", "[", i, "]", sep="")
}

Chill_Force$parameter<-new.names
mod.ranef2 <- full_join(mod.ranef2, Chill_Force)

Light_Force <- coef(modelhere2, prob=c(0.25, 0.75))$GEN.SPA[, c(1, 3:4), 7] %>%
  as.data.frame() %>%
  round(digits = 2) %>% 
  rename(mean = Estimate) %>%
  rename(`25%` = Q25) %>%
  rename(`75%` = Q75) %>%
  dplyr::select( mean, `25%`, `75%`) 
new.names<-NULL
for(i in 1:length(spp)){
  new.names[i]<-paste("Light:Force", "[", i, "]", sep="")
}

Light_Force$parameter<-new.names
mod.ranef2 <- full_join(mod.ranef2, Light_Force)

modoutput1 <- tidy(modelhere, prob=c(0.5))
modoutput2 <- tidy(modelhere2, prob=c(0.5))
