# Creating maps and assigning treatments for Budburst experiment 2016

# This script will assigns treatments and chamber locations for each cutting.
setwd("~/Documents/git/budchill")

library(gdata) # for read.xls
d <- read.xls("~/Documents/git/budchill/Field Sampling Dec 2015.xlsx")

allsp <- sort(unique(d[,"Species"]))

# Vib cas: only 4 inds sampled, still use

# Chilling treatment combinations
chill <- gl(4, 3, labels = c("chill1", "chill2","chill4","chill8"))
time <- gl(3, 1, length = 12, labels = c("time1", "time2","time3"))
treatcode <- paste(substr(chill, 6, 6), substr(time, 5, 5), sep = "_") 

chilltreat <- data.frame(chill, time, treatcode)


# Make twig-wise data frame. Assign treatments by individual, using chilltreat and nonchilltreat dataframes.

dx <- vector()

for(i in 1:nrow(d) ) { # i = 1


	xx <- paste(d[i,"Individual"], formatC(1:12, width = 2, flag = "0"), sep = "_")

	xx <- data.frame(xx)
	
	xx$sp <- substr(xx[,1], 1, 6)
	xx$rep <- substr(xx[,1], 7, 8)
	xx$ind <- substr(xx[,1], 1, 8)
	xx$twig <- substr(xx[,1], 13, 14)

	names(xx)[1] = "id"

	# Assign treatments. Randomize rows of treatment dataframes and apply to this individual
	xx <- data.frame(xx, chilltreat[sample(1:nrow(chilltreat)),])

	dx <- rbind(dx, xx)
	}

#write.csv(dx, paste("Budburst Twig Datasheet w Treatments ", Sys.Date(), ".csv", sep = ""), row.names=F)

# <<>> <<>> <<>> <<>> <<>> <<>> <<>> <<>> <<>> <<>> <<>> <<>> <<>> <<>> <<>> <<>> <<>> <<>> <<>> <<>>
# make maps. 
# <<>> <<>> <<>> <<>> <<>> <<>> <<>> <<>> <<>> <<>> <<>> <<>> <<>> <<>> <<>> <<>> <<>> <<>> <<>> <<>>

# Jan 2016 maps for twigs
# Divide by chill treatment.

jz <- vector()

for(i in unique(dx$chill)){ # i = "chill1" 
	jx <- dx[dx$chill == i,]
	# Randomize and make into sets of three; make sure to not have congeners in pairs.
	 
	# Positions: 52
	jx$position <- sample(rep(1:52, 3))
	
	jx <- jx[order(jx$position),]

	# Check to see if a pair has a conspecific. If so, repeat shuffling until no pairs

	jx2 <- vector() # to hold ones which do not have conspecifics in the beaker
	
	while(any(unlist(tapply(jx$sp, jx$position, duplicated))) & !(nrow(jx) <= 12 & nrow(jx) >= 3)){
		noconspecs <- !unlist(lapply(tapply(jx$sp, jx$position, duplicated), function(x) any(x)))

		jx2 <- rbind(jx2, jx[jx$position %in% names(noconspecs[noconspecs]),])
		jx <- jx[jx$position %in% names(noconspecs[!noconspecs]),]
		jx$position <- sample(jx$position)	
	
		jx <- jx[order(jx$position),]
				
		}
	if(nrow(jx) > 0){ jx <- rbind(jx2, jx) } else { jx <- jx2 }
	
	jx <- jx[order(jx$position),]
	
	jz <- rbind(jz, jx)
	}

# Write it out. Important: if the script is re-run, will get different assignments of treatments! So only do this once.
# Here using Sys.Date to prevent accidental overwriting of original csv.

# Adding in row and column values for each twig. Now with 52 places, need 8 rows and 9 columns. 
jz$row = rep(c(rep(sort(rep(1:8, 3)), 6), sort(rep(1:4, 3))), 4)
jz$col = rep( sort( c( rep(rep(1:7, 3), 4), rep(rep(1:6, 3), 4) ) ) , 4)

# jz[1:100, ]

write.csv(jz, paste("design/Budburst Twig Datasheet w Treatments ", Sys.Date(), ".csv", sep = ""), row.names=F)

# Make maps.

maplist <- list()

# 52 positions. 8 rows, 7 columns (last 4 cells of 7th col blank)

counter = 1
for(i in unique(jz$chill)){ # i = "chill1" 
	jx <- jz[jz$chill == i,]
	
	mat <- matrix(nrow = 8, ncol = 7, byrow = T)
	
	for(k in unique(jx$position)){ # k = 1
		mat[k] = paste(jx[jx$position == k,"id"], collapse = "\t")
		
		}
	
	maplist[[counter]] = mat
	
	counter = counter + 1
	}

write.csv(maplist, "design/Mapschill.csv")


##### Write out position-wise for labels

labvec <- vector()

for(i in unique(jz$chill)){ # i = "chill1"
	
	jx <- jz[jz$chill == i,]
	
	for(j in unique(jx$position)){ # j = 1
		xx <- as.character(jx[jx$position == j,"id"])
		if(length(xx) == 1) xx <- c(xx, " ", " ")
		if(length(xx) == 2) xx <- c(xx, " ")
		xx <- c(xx, as.character(jx[jx$position == j,"treatcode"][1]))
		
		labvec <- rbind(labvec, xx)
		}
		}

write.csv(labvec, "design/Avery Labels 3per.csv", row.names = F)

