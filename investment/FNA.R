### Started by Dan Feb 9 2021
### Part 1:### modeling FLS for FNA with measurement models in brms
#FLS ~ pdsi + flower traits + fruit traits
###Part 2: Focusing on prunocerasus with stan model
#based on 

rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
options(mc.cores = parallel::detectCores())


library(dplyr)
library("housingData")
library(stringr)
library("ncdf4")
library(raster)
library(ggrepel)
library(patchwork)
library(brms)
library("grid")
library(tidybayes)
library(bayesplot)
library("caper")
library(phytools)
library(geiger)


#-------------------------#
#------Part 1--------------#
#--------------------------#
setwd("~/Documents/git/proterant/investment")
load("FNA.Rda")
FNA<-read.csv("Data/cherry_data.csv") ## measurement and FLS data from FNA
##clean specicies name

FNA$species<-ifelse(FNA$species=="hortunlana","hortulana",FNA$species)
FNA$species<-ifelse(FNA$species=="gladulosa","glandulosa",FNA$species)
FNA$species<-ifelse(FNA$species=="speciosa","speciosa",FNA$species)
FNA$species<-ifelse(FNA$species=="fasiculata","fasciculata",FNA$species)


sort(unique(FNA$species))
#### get herbaria specimen coordinates
mid.herb<-read.csv("SymbOutput_2020-10-26_133412_DwC-A/occurrences.csv")
mid.herb<-filter(mid.herb,specificEpithet %in% unique(FNA$species)) ## filter to Species we have data 4
unique(mid.herb$specificEpithet)
## select rows that already have coordinate
pruno.ref<-filter(mid.herb,!is.na(decimalLatitude))
table(pruno.ref$stateProvince)
table(pruno.ref$specificEpithet)
pruno.ref$lon<-pruno.ref$decimalLongitude
pruno.ref$lat<-pruno.ref$decimalLatitude

##now add county level coordiantes
#### two without coordinates
pruno.unref<-filter(mid.herb,is.na(decimalLatitude)) ## filter entries with no lat/lon

### prep county data
geoCounty$rMapState<-str_to_title(geoCounty$rMapState) ### centriod coordinates for all US counties
colnames(geoCounty)[6]<-"stateProvince" 

colnames(geoCounty)[c(2,7)]<-c("County_name","county") ## MATCH NAMES

pruno.unref<-filter(pruno.unref,id!=17453277) #remvoe problematic canadian entry
pruno.unref$county<-tolower(pruno.unref$county) ## make ours lower case



pruno.unref$county<-gsub("county","",pruno.unref$county) ### get rid of extenious coounty info
pruno.unref$county<-gsub("()","",pruno.unref$county) # ""
pruno.unref$county<-gsub("co.","",pruno.unref$county,fixed = TRUE) # ""




pruno.unref<-dplyr::left_join(pruno.unref,geoCounty,by=c("county","stateProvince"))
pruno.unref<-dplyr::select(pruno.unref,-County_name,-fips, -state)
intersect(colnames(pruno.unref),colnames(pruno.ref))

pruno.ref$geoclass<-"geo_referenced"
pruno.unref$geoclass<-"county_centroid"
pruno.unref<-filter(pruno.unref,!is.na(lon)) # select only entries with lat lon
pruneo<-rbind(pruno.ref,pruno.unref)


pruneo<-filter(pruneo,lat>0)
pruneo<-filter(pruneo,lon<(-50))
ggplot(pruneo,aes(lon,lat))+geom_point(aes(color=specificEpithet))

##now add pdsi data
palmer.b <- brick("Data/lbda-v2_kddm_pmdi_2017.nc")

lonpoints<-pruneo$lon # make vector of prunus coordinates
latpoints<-pruneo$lat #
extract.pts <- cbind(lonpoints,latpoints)

mean.prunus <-calc(palmer.b, fun = mean,na.rm=TRUE) #average palmer drought index acrosss time
sd.prunus <-calc(palmer.b, fun = sd,na.rm=TRUE) #
#min.prunus <-calc(palmer.b, fun = min,na.rm=TRUE) 
ext<-raster::extract(mean.prunus,extract.pts,method="simple")
ext2<-raster::extract(sd.prunus,extract.pts,method="simple")
ext3<-raster::extract(min.prunus,extract.pts,method="simple")
pruneo$pdsi<-ext
pruneo$pdsi.sd<-ext2
pruneo$pdsi.min<-ext3


pdsi<- pruneo %>% dplyr::group_by(specificEpithet) %>% dplyr::summarise(meanpdsi=mean(pdsi,na.rm=TRUE),sdmean=sd(pdsi,na.rm=TRUE),minpdsi=mean(pdsi.min,na.rm=TRUE),sdmin=sd(pdsi.min,na.rm=TRUE))
colnames(pdsi)[1]<-"species"
setdiff(FNA$species,pdsi$species) ## lose speciosa


## combine data
FNA<-left_join(FNA,pdsi)
##quick plots
cor(FNA$petal_low,FNA$petal_high) # .88
cor(FNA$inflor_low,FNA$inflor_high) #.86 can probably use jsut one or the other in the
cor(FNA$inflor_low,FNA$petal_low) #-.35 right direction but not as correlated as trade off would suggest
cor(FNA$inflor_high,FNA$petal_high) #-.31
cor(FNA$fruit_low,FNA$fruit_high) #.88
##Make FLS numeric for oridnal modeling
unique(FNA$FLS)
FNA$FLSnum<-NA
FNA$FLSnum[FNA$FLS=="before"]<-1
FNA$FLSnum[FNA$FLS=="before/with"]<-2
FNA$FLSnum[FNA$FLS=="with"]<-3
FNA$FLSnum[FNA$FLS=="after"]<-4
FNA$FLSnum<-as.integer(FNA$FLSnum)
##quick plots
class(FNA$FLSnum)
plot(FNA$FLSnum,FNA$inflor_low)

##basic model
# FLSnum~ inflorlow*petal_low+me(pdsi,sd)+fruit_low +phylo someday

#rough sd and mean calculation for other predictons
#FNA$petalmean<-(FNA$petal_low+FNA$petal_high)/2
#FNA$petalsd<-(FNA$petal_high-FNA$petal_low)/4 #"rule of ranges"


#FNA$inflormean<-(FNA$inflor_low+FNA$inflor_high)/2
#FNA$inflorsd<-(FNA$inflor_high-FNA$inflor_low)/4 #"rule of ranges"

#FNA$fruitmean<-(FNA$fruit_low+FNA$fruit_high)/2
#FNA$fruitsd<-(FNA$fruit_high-FNA$fruit_low)/4 #"rule of ranges"

zscore <- function(x){(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)}
#FNA$fruit.z.mean<-range01(FNA$fruitmean)
#FNA$fruit.z.sd<-range01(FNA$fruitsd)

#FNA$inflor.z.mean<-range01(FNA$inflormean)
#FNA$inflor.z.sd<-range01(FNA$inflorsd)

#FNA$petal.z.mean<-range01(FNA$petalmean)
#FNA$petal.z.sd<-range01(FNA$petalsd)

#zscore all predictor
FNA$meanpdsi.z<-zscore(FNA$meanpdsi)
FNA$minpdsi.z<-zscore(FNA$minpdsi)
FNA$petal.z<-zscore(FNA$petal_high)
FNA$inflor.z<-zscore(FNA$inflor_high)
FNA$fruit.z<-zscore(FNA$fruit_high)
#FNA$cold.z<-zscore(FNA$meanT)

#########################################################################
#######Gaussian models with and without phylogenies#####################
#########################################Note: they dont converge with phylos
#check colinearity
cor(FNA$petal.z,FNA$inflor.z)
cor(FNA$fruit.z,FNA$inflor.z)
cor(FNA$fruit.z,FNA$meanpdsi.z,use="complete.obs")
cor(FNA$petal.z,FNA$meanpdsi.z,use="complete.obs")
cor(FNA$inflor.z,FNA$meanpdsi.z,use="complete.obs")

FNA$logFLS<-log(FNA$FLSnum)

FNAgaus<-brm(logFLS~petal_high*inflor_high+fruit_high+meanpdsi,
               data=FNA,warmup=3000,iter=4000)

FNAgaus.z<-brm(logFLS~petal.z*inflor.z+fruit.z+meanpdsi.z,
             data=FNA,warmup=3000,iter=4000)

cherrytree<-read.tree("pruned_4_modeling.tre")
is.ultrametric(cherrytree)
FNA.small<-filter(FNA,species %in% cherrytree$tip.label)
A <- ape::vcv.phylo(cherrytree)

FNAgaus.z.small<-brm(logFLS~petal.z*inflor.z+fruit.z+meanpdsi.z,
               data=FNA.small,warmup=3000,iter=4000)

get_prior(logFLS~meanpdsi.z+petal.z+fruit.z +inflor.z+ (1|gr(species, cov = A)),
          data=FNA.small,
          data2 = list(A = A))

priorz<-c(prior_string("student_t(3,1,4)",class="b"),
          prior_string("student_t(3,0,3)",class="Intercept"),
          prior_string("student_t(3,0,2.5)",class="sd"),
          prior_string("student_t(3,0,3)",class="sigma"))
plot(FNA.small$fruit.z)

FNAgaus.phylo<-brm(logFLS~meanpdsi+petal_high+fruit_high+ (1|gr(species, cov = A)),
                   data=FNA.small,
                   data2 = list(A = A),control=list(adapt_delta=0.99),
                   warmup=5000,iter=6000)
fixef(FNAgaus.phylo,probs=c(.025,0.25,.75,.975))
tidybayes::get_variables(FNAgaus.phylo)
goo<-spread_draws(FNAgaus.phylo,b_meanpdsi, b_petal_high,b_fruit_high[condition,term])
head(goo)
FNAgaus.phylo2<-brm(logFLS~inflor_high+fruit_high+meanpdsi+ (1|gr(species, cov = A)),
                   data=FNA.small,
                   data2 = list(A = A),control=list(adapt_delta=0.99),
                   warmup=3000,iter=4000)

launch_shinystan(FNAgaus.phylo)
fixef(FNAgaus.z)
fixef(FNAgaus.z,probs = c(0.25,.75))
fixef(FNAgaus.z.small,probs = c(0.25,.75))
fixef(FNAgaus.phylo,probs = c(0.25,.75))
fixef(FNAgaus.phylo2,probs = c(0.25,.75))
summary()

fixef(FNAgaus.phylo,probs = c(0.25,.75))

priorz<-get_prior(FLSnum~petal.z+inflor.z+meanpdsi.z + (1|gr(species, cov = A)),
          data=FNA.small,
          data2 = list(A = A))

priorz$prior[5]<-set_prior("student_t(3, 3, 1.5)")

FNAgaus.phylo<-brm(FLSnum~meanpdsi.z + (1|gr(species, cov = A)),
                   data=FNA.small,
                   data2 = list(A = A),
                   warmup=6500,iter=8000,control=list(adapt_delta=.995))

pp_check(FNAgaus.phylo)

summary(FNAgaus.phylo)
### try it in phylolm
mytree.names<-cherrytree$tip.label


###format the data in the same order as the tree
final.df<-FNA.small[match(mytree.names, FNA.small$species),]
namelist2<-final.df$species
namelist2==mytree.names
final.df$species== mytree.names

final.df<-  final.df %>% remove_rownames() %>% column_to_rownames(var="species")
library(tibble)
library("phylolm")


FNAgaus.phylo.lm<-phylolm(logFLS~petal.z*inflor.z+fruit.z+meanpdsi.z,data=final.df,phy=cherrytree,boot = 9999)

final.df2<-FNA.small[match(mytree.names, FNA.small$species),]
FNAgaus.phylo.brm<-brm(logFLS~petal.z*inflor.z+fruit.z+meanpdsi.z+(1|gr(species, cov = A)),
                    data=final.df2,
                    data2 = list(A = A),control=list(adapt_delta=0.99),
                    warmup=3000,iter=4000)


summary(FNAgaus.phylo.lm)
simplelm<-lm(logFLS~petal.z+inflor.z+fruit.z+meanpdsi.z,data=final.df)
summary(simplelm)

summary(FNAgaus.phylo.lm)
summary(simplelm)
fixef(FNAgaus.z)



####################################################################################################
#########################################################################
#######Oridinal models with and without phylogenies#####################
#########################################Note:

##########
FNAordz<-brm(FLSnum~petal.z*inflor.z+fruit.z+meanpdsi.z,
                   data=FNA,
                   family=cumulative("logit"),warmup=3000,iter=4000)


get_variables(FNAordz)
output<-FNAordz %>%
   spread_draws(b_petal.z,b_inflor.z ,b_fruit.z,b_meanpdsi.z,`b_petal.z:inflor.z`)
colnames(output)
output <-output %>% tidyr::gather("var","estimate",4:8)


FNAordz.phylo<-brm(FLSnum~petal.z+fruit.z+meanpdsi.z + (1|gr(species, cov = A)),
             data=FNA.small,
             family=cumulative("logit"),
             data2 = list(A = A),control=list(adapt_delta=0.9),
             warmup=5000,iter=6000)
fixef(FNAordz.phylo,probs = c(.025,.25,.75,.975))
FNA.small$specificEpithet<-FNA.small$species

FNA.small$meanpdsi.z<-zscore(FNA.small$meanpdsi)
FNA.small$inflor.z<-zscore(FNA.small$inflor_high)
FNA.small$fruit.z<-zscore(FNA.small$fruit_high)







fixef(FNAordz.phylo,probs = c(.25,.75))
fixef(FNAordz.phylo2,probs = c(.05,.25,.75,.95))
conditional_effects(FNAordz.phylo2,prob = .5,categorical=F)
conditional_effects(FNAgaus.phylo,prob = .5)


###reverseFLS number so on same scale as Plums
FNA.small$FLSnum[FNA.small$FLS=="before"]<-4
FNA.small$FLSnum[FNA.small$FLS=="before/with"]<-3
FNA.small$FLSnum[FNA.small$FLS=="with"]<-2
FNA.small$FLSnum[FNA.small$FLS=="after"]<-1
FNA.small$FLSnum<-as.integer(FNA.small$FLSnum)


###this is the main model
FNAordz.phylo2<-brm(FLSnum~inflor.z*meanpdsi.z +(1|specificEpithet)+(1|gr(species, cov = A)),
                    data=FNA.small,
                    family=cumulative("logit"),
                    data2 = list(A = A),control=list(adapt_delta=0.99),
                    warmup=6000,iter=8000)




FNAordz.small<-brm(FLSnum~petal.z*inflor.z+fruit.z+meanpdsi.z,
             data=FNA.small,
             family=cumulative("logit"),warmup=3000,iter=4000)

FNAgaus.small<-brm(FLSnum~petal.z*inflor.z+fruit.z+meanpdsi.z,
                   data=FNA.small,warmup=3000,iter=4000)
#HERE
fixef(FNAordz.phylo, probs = c(.25,0.75))
fixef(FNAordz.small,probs = c(.25,0.75))
fixef(FNAordz,probs = c(.25,0.75))


fixef(FNAgaus.phylo)
fixef(FNAgaus.small)
fixef(FNAordz)


round(fixef(FNAordz.phylo2, probs = c(.055,.25,0.75,.945)),2)
posterior <- as.matrix(FNAordz)
colnames(posterior)
mcmc_areas(posterior,
           pars = c("b_petal.z","b_inflor.z"  ,"b_fruit.z","b_meanpdsi.z","b_petal.z:inflor.z"),
           prob = 0.5)


get_variables(FNAordz.phylo2)
output<-FNAordz.phylo2 %>%
   spread_draws(b_inflor.z ,b_meanpdsi.z,`b_inflor.z:meanpdsi.z`)
colnames(output)
output <-output %>% tidyr::gather("var","estimate",4:6)



#jpeg("Plots/fullprunus_mus.jpeg",width=9, height=7, units = "in",res=300)
pottymu<-ggplot(output,aes(y = var, x = estimate)) +
   stat_pointinterval(.width=c(.5,.89))+ggthemes::theme_few()+
   geom_vline(xintercept=0,linetype="dashed")+xlim(-50,10)+
   scale_y_discrete(limits = c("b_inflor.z:meanpdsi.z","b_meanpdsi.z","b_inflor.z"),labels=c("PDSI X infloresence size","PDSI","inflorescence size"))+
   ylab("")+xlab("standardized effect size estimate")

#dev.off()
#FNAord.traits<-brm(FLSnum~me(petalmean,petalsd)+me(inflormean,inflorsd)+
 #                    me(fruitmean,fruitsd),
  #                 data=FNA,
   #                family=cumulative("logit"),warmup=3000,iter=4000,control = list(adapt_delta=.95))

launch_shinystan(FNAordz.error)


#get_prior(FLSnum~petal.z.mean+inflor.z.mean+
 # me(meanpdsi.z,sdpdsi.z)+fruit.z.mean,
#data=FNA)

#bprior <- c(prior(normal(0,1), class = b),
 #           prior(normal(0,1), class = meanme ),
  #          prior(uniform(0,1), class = sdme), 
   #               prior(normal(0,1), class = Intercept))

#FNAordz.error1<-brm(FLSnum~petal.z.mean+inflor.z.mean+
 #                    meanpdsi.z+me(fruit.z.mean,fruit.z.sd),
  #                 data=FNA ,save_mevars = TRUE,
   #                family=cumulative("logit"),warmup=3000,iter=4000 )

#launch_shinystan(FNAordz.error1)
p1<-plot(conditional_effects(FNAordz.phylo2, "meanpdsi.z", ordinal = TRUE,prob = .5,plot=FALSE))
p3<-plot(conditional_effects(FNAordz.phylo2, "inflor.z", ordinal = TRUE,prob=.5,plot=FALSE))

conditions <- make_conditions(FNAordz.phylo2, "inflor.z")
p4<-plot(conditional_effects(FNAordz.phylo2, "meanpdsi.z",conditions=conditions,ordinal = TRUE,prob=.5,plot=FALSE))
range(FNA.small$inflor.z)
p1<-p1[[1]]+ggthemes::theme_few()+scale_y_discrete(name="FLS",labels=c("flowers after leaves","flowers with leaves","flowers before/with leaves","flowers before leaves"))+xlab("PDSI")

#p2<-p2[[1]]+ggthemes::theme_few()+scale_color_manual(values=c("hotpink","orange","lightgreen","darkgreen"))+scale_fill_manual(values=c("hotpink","orange","lightgreen","darkgreen"))
p3<-p3[[1]]+ggthemes::theme_few()+ylab("")+xlab("inflorescence size")+theme(axis.text.y=element_blank(),axis.ticks.y = element_blank())

p4<-p4[[1]]+ggthemes::theme_few()+xlab("inflorescence size")+scale_y_discrete(name="FLS",labels=c("flowers before leaves","flowers before/with leaves","flowers with leaves","flowers after leaves"))

+scale_color_manual(name="FLS",labels=c("flowers before leaves","flowers before/with leaves","flowers with leaves","flowers after leaves"),values=c("hotpink","orange","lightgreen","darkgreen"))+
   scale_fill_manual(name="FLS",labels=c("flowers before leaves","flowers before/with leaves","flowers with leaves","flowers after leaves"),values=c("hotpink","orange","lightgreen","darkgreen"))+xlab("fruit size")

potty<-ggpubr::ggarrange(p1,p3,common.legend = TRUE,ncol=2,legend="bottom",widths = c(.8,.5))

jpeg("Plots/whatreviwerswant/fullprunus_4manu.jpeg",width=9, height=8, units = "in",res=300)
ggpubr::ggarrange(pottymu,potty,nrow=2,labels=c("a)","b)"))
dev.off()
#### run this just on prunocerasus
pruno<-c("alleghaniensis","angustifolia","americana" ,"gracilis","geniculata","hortulana" ,"maritima",
         "mexicana","murrayana","munsoniana","nigra","rivularis","umbellata","subcordata","texana" )


FNA.pruno<-dplyr::filter(FNA, species %in% pruno)

FNAordz.pruno<-brm(FLSnum~petal.z+fruit.z+cold.z+meanpdsi.z,
             data=FNA.pruno,
             family=cumulative("logit"),warmup=3000,iter=4000)

fixef(FNAordz.pruno)


p11<-plot(conditional_effects(FNAordz.pruno, "meanpdsi.z", categorical = TRUE,probs = c(.25,.75),plot=FALSE))
p55<-plot(conditional_effects(FNAordz.pruno, "cold.z", categorical = TRUE,probs = c(.25,.75),plot=FALSE))
p22<-plot(conditional_effects(FNAordz.pruno, "petal.z", categorical = TRUE,probs = c(.25,.75),plot=FALSE))
p44<-plot(conditional_effects(FNAordz.pruno, "fruit.z", categorical = TRUE,probs = c(.25,.75),plot=FALSE))

p11<-p11[[1]]+ggthemes::theme_few()+scale_color_manual(values=c("hotpink","orange","lightgreen","darkgreen"))+scale_fill_manual(values=c("hotpink","orange","lightgreen","darkgreen"))
p22<-p22[[1]]+ggthemes::theme_few()+scale_color_manual(values=c("hotpink","orange","lightgreen","darkgreen"))+scale_fill_manual(values=c("hotpink","orange","lightgreen","darkgreen"))
p44<-p44[[1]]+ggthemes::theme_few()+scale_color_manual(values=c("hotpink","orange","lightgreen","darkgreen"))+scale_fill_manual(values=c("hotpink","orange","lightgreen","darkgreen"))
p55<-p55[[1]]+ggthemes::theme_few()+scale_color_manual(values=c("hotpink","orange","lightgreen","darkgreen"))+scale_fill_manual(values=c("hotpink","orange","lightgreen","darkgreen"))

jpeg("Plots/FNA_mean_prunocerasus.jpeg")
ggpubr::ggarrange(p11,p55,p22,p44,common.legend = TRUE)
dev.off()
fixef(FNAordz.pruno,probs = c(.25,.75))

cor(FNA$petal_low,FNA$meanpdsi,use = "complete.obs")
cor(FNA$minpdsi,FNA$meanpdsi,use = "complete.obs")

cor(FNA.pruno$minpdsi,FNA.pruno$meanpdsi,use = "complete.obs")
cor(FNA.pruno$meanT,FNA.pruno$meanpdsi,use = "complete.obs")





#treee<-read.tree("..//input/Vascular_Plants_rooted.dated.tre")
#names.intree<-treee$tip.label
#namelist<-FNA$name

#which(names.intree%in%namelist)
#length(which(namelist%in%names.intree))
#length(which(!namelist%in%names.intree))

save.image("FNA.Rda")


