###This should be sourced in concept.paper.figures to make a more restricted dataset to test flowering. Made Jan 2018
median(mich.data$flo_time)
mean(mich.data$flo_time)
earlymich<-filter(mich.data,flo_time<=5.5)

###prune tree
mich.data<-rownames_to_column(mich.data, "name")
namelist<-earlymich$name
names.intree<-mich.tree$tip.label

##Prune the tree
to.prune<-which(!names.intree%in%namelist)
earlytree<-drop.tip(mich.tree,to.prune)

earlytree$tip.label==earlymich$name

###For some reason I need to runthis twice when sourcing
earlymich<-filter(mich.data,flo_time<=5.5)

###prune tree
#mich.data<-rownames_to_column(mich.data, "name")
namelist<-earlymich$name
names.intree<-mich.tree$tip.label

##Prune the tree
to.prune<-which(!names.intree%in%namelist)
earlytree<-drop.tip(mich.tree,to.prune)

earlytree$tip.label==earlymich$name

#recenter
earlymich$height_cent<-(earlymich$heigh_height-mean(earlymich$heigh_height))/(2*sd(earlymich$heigh_height))
earlymich$fruit_cent<-(earlymich$fruiting-mean(earlymich$fruiting))/(2*sd(earlymich$fruiting))
earlymich$flo_cent<-(earlymich$flo_time-mean(earlymich$flo_time))/(2*sd(earlymich$flo_time))
earlymich$pol_cent<-(earlymich$pol-mean(earlymich$pol))/(2*sd(earlymich$pol))
earlymich$av_fruit_time_cent<-(earlymich$av_fruit_time-mean(earlymich$av_fruit_time))/(2*sd(earlymich$av_fruit_time))

earlymich<-  earlymich %>% remove_rownames %>% column_to_rownames(var="name")

earlycent<-phyloglm(pro2~pol+height_cent+flo_cent+fruit_cent+shade_bin,earlymich, earlytree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                    start.beta=NULL, start.alpha=NULL,
                    boot=50,full.matrix = TRUE)
summary(earlycent)
