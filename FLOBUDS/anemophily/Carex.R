rm(list=ls())
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
#rstan_options(auto_write = TRUE)
graphics.off()
library(dplyr)
library(brms)

setwd("~/Desktop/exploratory/")

original<-read.csv("gamma_grass_original.csv")
checky<-as.data.frame(colnames(original))
original<-dplyr::select(original,stateProvince,gbifID,references,	 verbatimElevation)
if(FALSE){
d<-read.csv("gamma_grass.csv")
colnames(d )



d<-dplyr::select(d,species,stateProvince,eventDate,year,month,day,decimalLatitude,decimalLongitude,references)

tab<-d %>% group_by(stateProvince) %>% count()%>% arrange(-n)

TabeE<-filter(tab, stateProvince %in% c("Florida","Georgia","South Carolina","North Carolina","Virginia","Maryland"))


TabeC<-filter(tab, stateProvince %in% c("Texas","Oklahoma","Kansas","Nebraska"))
d$PIST<-NA
d$STAM<-NA
d <- d %>%
  dplyr::arrange(match(stateProvince, tab$stateProvince))

write.csv(d,"gamma_grass.csv")
}
d<-read.csv("gamma_grass_working.csv")

#d<-read.csv("carex_data - carex_data.csv")
#unique(d$specificEpithet)
#d<-dplyr::filter(d,specificEpithet %in% c("plantaginea","platyphylla"))

d<-filter(d, stateProvince=="North Carolina")
if(FALSE){
#d<-dplyr::filter(d,flowers_visable==1)
d<-dplyr::filter(d, d$PIST==1 & d$STAM==0 | d$PIST==1 & d$STAM==1)


d$response<-ifelse(d$PIST==1 & d$STAM==0,1,0)

table(d$response)

d <- d %>%
  mutate(
    date = make_date(year, month, day),
    doy  = yday(date)
  )

d$year2<-ifelse(d$year<1980,1980,d$year)

fit1<-brm(response~doy, family=bernoulli, data=d)

fit2<-brm(response~decimalLatitude+doy*year2, family=bernoulli, data=d)

brms::conditional_effects(fit1)

}

d <- d %>%
  mutate(
    date = make_date(year, month, day),
    doy  = yday(date)
  )

colnames(d)
d<-left_join(d,original)


d<-gather(d, "phase", "response",10:11)

d$response<-ifelse(d$response==1.5,0,d$response)
d$response<-ifelse(d$response==0.5,0,d$response)

unique(d$phase)


library(ggplot2)


goo<-filter(d,response==1)


ggplot(goo,aes(doy))+geom_density(aes(fill=phase),position="dodge",alpha=0.5)

ggplot(goo, aes(elevation,doy))+geom_smooth(method="lm",aes(color=phase))

library(prism)
library(terra)   
library(lubridate)
library(purrr)

prism_set_dl_dir("~/prism_data")

years <- sort(unique(goo$year))

df<-goo

walk(years, function(y) {
  
  message("Downloading PRISM for year ", y)
  
  max_doy <- max(df$doy[df$year == y], na.rm = TRUE)
  
  if (!is.finite(max_doy) || max_doy < 1) {
    warning("Skipping year ", y, ": invalid DOY")
    return(NULL)
  }
  
  start_date <- as.Date(sprintf("%d-01-01", y))
  end_date   <- start_date + max_doy - 1
  
  get_prism_dailys(
    type = c("tmin", "tmax"),
    minDate = start_date,
    maxDate = end_date,
    keepZip = FALSE
  )
})





mod1<-brm(doy~decimalLatitude+phase, data=goo)
conditional_effects()

goober$year2<-ifelse(goober$year<1980,1980,goober$year)
library(brms)

table(goo$phase,goo$year2)

prior <- c(
  prior(normal(0, 50), class = "b"),      # stronger regularization
  prior(normal(365, 100), class = "Intercept")
)



mod2<-brm(doy~decimalLatitude+(decimalLatitude|phase), data=goober)

coef(mod2
     )
conditional_effects(mod2)



