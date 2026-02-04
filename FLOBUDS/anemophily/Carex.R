rm(list=ls())
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
#rstan_options(auto_write = TRUE)
graphics.off()
library(dplyr)
library(brms)
library(lubridate)

setwd("~/Documents/git/proterant/FLOBUDS/anemophily/")

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
d<-read.csv("gamma_grass - gamma_grass.csv")

#d<-read.csv("carex_data - carex_data.csv")
#unique(d$specificEpithet)
#d<-dplyr::filter(d,specificEpithet %in% c("plantaginea","platyphylla"))

d<-filter(d, stateProvince %in% c("North Carolina","Virginia","Maryland","Illinois","Georgia","New York", "Pennsylvania","South Carolina"))
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

#### Hearns package
install.packages("copula")
devtools::install_github("david-j-hearn/phenoCollectR", dependencies = TRUE)

remotes::install_github("david-j-hearn/phenoCollectR")
library(phenoCollectR)
installed.packages()[, "Package"]

ls("package:phenoCollectR")




library(phenoCollectR)
vars<-c("phase","decimalLatitude")
data <- preparePhenologyData(dataFile=goo, responseVariableName="doy", onsetCovariateNames=vars, durationCovariateNames=vars, removeOutliers=TRUE)




ggplot(goo,aes(doy))+geom_density(aes(fill=phase),position="dodge",alpha=0.5)

ggplot(goo, aes(decimalLatitude,doy))+geom_smooth(method="lm",aes(color=phase))


library(daymetr)
library(prism)
library(terra) 
library(dplyr)
library(lubridate)
library(purrr)


range(d$year)

##get monthy temperature atara
prism_set_dl_dir("~/prism_data")

years <- sort(unique(goo$year))

get_prism_monthlys(
  type = "tmean",
  years = years,
  keepZip = FALSE
)


dirs <- prism_archive_ls()

# Keep monthly tmean datasets
dirs <- dirs[grepl("tmean", dirs)]

# Build full path to the .bil files
bil_files <- file.path(
  prism_get_dl_dir(),
  dirs,
  paste0(dirs, ".bil")
)

# Sanity check (important)
bil_files <- bil_files[file.exists(bil_files)]





# Load rasters
tmean <- rast(bil_files)




if(FALSE){

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
}

# Add time index
dates <- as.Date(
  paste0(
    sub(".*_(\\d{4})(\\d{2})$", "\\1", dirs),
    "-",
    sub(".*_(\\d{4})(\\d{2})$", "\\2", dirs),
    "-15"
  )
)
time(tmean) <- dates


#dates <- as.Date(paste0(tmean_files$year, "-", tmean_files$month, "-15"))
time(tmean) <- dates


tmean_annual <- tapp(tmean, format(time(tmean), "%Y"), mean, na.rm = TRUE)


df <- dplyr::select(goo,decimalLatitude,decimalLongitude, year)

pts <- vect(df, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
pts <- project(pts, crs(tmean_annual))

vals <- extract(tmean_annual, pts)


names(tmean_annual) <- gsub("^X", "", names(tmean_annual))
years_available <- as.integer(names(tmean_annual))

df_out <- df |>
  mutate(
    col_id = match(year, years_available),
    tmean = vals[cbind(seq_len(nrow(vals)), col_id)]
  ) |>
  select(decimalLongitude, decimalLatitude, year, tmean)


goober<-goo
goober$meanT<-df_out$tmean

ggplot(goober, aes(meanT,doy))+geom_smooth(method="lm",aes(color=phase))


prior <- c(
  prior(normal(0, 100), class = "b"),      # stronger regularization
  prior(normal(200, 100), class = "Intercept")
)
goo$lat_cent<-goo$decimalLatitude-mean(goo$decimalLatitude)
goo$year_cent<-goo$year-mean(goo$year)

mod1<-brm(doy~lat_cent+phase+year_cent+lat_cent:phase+year_cent:phase,control=list(adapt_delta=0.99),prior = prior,warmup=3000, iter=4000,data=goo)
conditional_effects(mod1,prob = .5)

ggplot(goober,aes(meanT))+geom_histogram()
ggplot(goo,aes(year))+geom_histogram()
ggplot(goo,aes(lat_cent))+geom_histogram()
 
goober<-filter(goober,meanT<30)      
mod2<-brm(doy~meanT*phase,control=list(adapt_delta=0.99),prior=prior,warmup=3000, iter=4000,data=goober)
conditional_effects(mod2,prob = .5)

get_prior(mod2)

priorz <- c(
  prior(normal(0, 100), class = "b"),      # stronger regularization
  prior(normal(200, 100), class = "Intercept"),
  prior(normal(0,100), class = "sd"),
  prior(normal(0,100), class = "sigma")
)
conditional_effects(mod2,prob = .5)

fixef(mod2,prob=.5)





get_daymet <- function(lat, lon, year) {
  download_daymet(
    
    lat = lat,
    lon = lon,
    start = year,
    end = year,
    internal = TRUE,
    silent = TRUE
  )$data
}

?pmap()
site_years <- goober %>%
  distinct(decimalLatitude, decimalLongitude, year)

site_years<-site_years[c(1:212,214:523),]
site_years<-site_years[c(1:212,214:348,351:523),]

daymet_data <- site_years %>%
  mutate(
    daymet = purrr::pmap(
      list(decimalLatitude, decimalLongitude, year),
      get_daymet
    )
  ) %>%
  unnest(daymet)


Tbase <- 10

daymet_data <- daymet_data %>%
  mutate(
    gdd_daily = pmax(
      ((tmax..deg.c. + tmin..deg.c.) / 2) - Tbase,
      0
    )
  ) %>%
  arrange(lat, lon, year, yday) %>%
  group_by(lat, lon, year) %>%
  mutate(gdd_cum = cumsum(gdd_daily)) %>%
  ungroup()

gdd_events <- events %>%
  left_join(
    daymet_data %>%
      select(lat, lon, year, yday, gdd_cum),
    by = c("lat", "lon", "year", "event_doy" = "yday")
  )








goober$year2<-ifelse(goober$year<1980,1980,goober$year)
library(brms)

table(goo$phase,goo$year2)





mod2<-brm(doy~decimalLatitude+(decimalLatitude|phase), data=goober)

coef(mod2
     )
conditional_effects(mod2)



