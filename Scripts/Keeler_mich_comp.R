###run Keeler.code. then line 24-26 from hyst_analysis
anthyK<-filter(final.df, name %in% namelistK)

##Prune the tree
names.intree<-pruned.by.anthy$tip.label

# list of my species myspecies
namelist<-unique(anthyK$name)

##Prune the tree
to.prune<-which(!names.intree%in%namelist)
pruned.by.anthyK<-drop.tip(pruned.by.anthy,to.prune)
pruned.by.anthyK$tip.label

###models
anthyK<-  anthyK %>% remove_rownames %>% column_to_rownames(var="name")

full.mod<-glm(pro~pol+heigh_height+shade_bin+av_fruit_time+flo_type+flo_time,family = binomial(link="logit"),data=anthyK)
summary(full.mod)

full.modB<-phyloglm(pro~pol+heigh_height+shade_bin+av_fruit_time+flo_time ,anthyK, pruned.by.anthyK, method = "logistic_MPLE", btol = 40, log.alpha.bound =4,
                    start.beta=NULL, start.alpha=NULL,
                    boot=10,full.matrix = TRUE)
summary(full.modB)

##### Okay, pollination syndrome drops out which is a good sign for analysis


## Now trying to figure out what av_fruit_time is doing
cor(anthyK$av_fruit_time, anthyK$flo_time)
plot(anthyK$av_fruit_time, anthyK$flo_time)

###center these variable
anthyK$centered_flo<-anthyK$flo_time/mean(anthyK$flo_time)
anthyK$centered_fruit<-anthyK$av_fruit_time/mean(anthyK$av_fruit_time)
ggplot(anthyK, aes(centered_flo,centered_fruit))+geom_point(aes(col=pro))


#### Do it for main dataset
final.df$centered_flo<-final.df$flo_time/mean(final.df$flo_time)
final.df$centered_fruit<-final.df$av_fruit_time/mean(final.df$av_fruit_time)

cor(final.df$centered_flo,final.df$centered_fruit)

ggplot(final.df, aes(centered_flo,centered_fruit))+geom_point(aes(col=pro))
ggplot(final.df, aes(flo_time,av_fruit_time))+geom_point(aes(col=pro))
ggplot(final.df, aes(flo_time,av_fruit_time))+geom_point(aes(col=Family))

#remove oak outlyers
outlyer<- filter(final.df, av_fruit_time<13)
cor(outlyer$flo_time,outlyer$av_fruit_time)

protan<-filter(final.df, pro==1)
cor(protan$flo_time,protan$av_fruit_time)
cor(protan$centered_flo,protan$centered_fruit)





