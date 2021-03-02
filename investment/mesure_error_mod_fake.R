rm(list=ls()) 
options(stringsAsFactors = FALSE)


N <- 1000
dat <- data.frame(
  y = rnorm(N), sdy = abs(rnorm(N, 1)), 
  x = rnorm(N), sdx = abs(rnorm(N, 1)))

#fit a simple error-in-variables model 


fit1 <- brm(y|mi(sdy) ~ me(x, sdx), data = dat, 
            save_mevars = TRUE) #runss but baby

fit2<-brm(y~x,data = dat)

fixef(fit2)
fixef(fit1)

load("Input/PrunusFLSs.Rda")

pdsi.group<-data.frame(species=mod.pdsi.out$species, pdsi=mod.pdsi.out$Estimate,pdsi.error=mod.pdsi.out$SE)
fls.group<-data.frame(species=mod2a.scale.out$species, fls=mod2a.scale.out$Estimate,fls.error=mod2a.scale.out$SE)

group<-left_join(fls.group,pdsi.group)

fit.fls <- brm(pdsi|mi(pdsi.error) ~ me(fls, fls.error), data = group, 
            save_mevars = TRUE)
prior_summary(fit.fls)

fixef(fit.fls,probs = c(.025,.25,.75,.975))
pp_check(fit.fls,nsamples = 500)
pp_check(fit.fls.noerr,nsamples = 500)

fit.fls.noerr <- brm(pdsi ~ fls, data = group, 
               save_mevars = TRUE)



fixef(fit.fls.noerr,probs = c(.025,.25,.75,.975))
