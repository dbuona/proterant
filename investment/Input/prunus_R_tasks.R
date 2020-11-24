## house keeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()

# Task 1: a) Set working directory
#         b) read in .csv file and call it "d"

#Task 2: explore the data
#     a) use functions nrow(), dim() colnames() to explore the data
      #b) how would you change the columnname of column 3 from "catalogNumber" to "barcode"...dooooo it.

#Task 3: What percetange of the datasheet has "Y" for d$flowers?
#         What percentage of each species has "Y" for d$species 
#     hint: use function table() for each


#Task 4: add a column (I'll help you with this one), add a column for "fruits"
# first explor the ifelse() function 
?ifelse()

d$fruits<- ifelse(d$bbch.f %in% c(89,NA),"Y","N") ## be prepare to explain to me what this funciton is doing here
d$flowers<-ifelse(d$fruits,"N",d$flowers)

#Task 5: plotting with ggplot2
#   a) load ggplot2 library
# basic ggplot syntax: ggplot(data.name,aes(x variable, y variable))+
#geom_object() A geom object is where you specify what kind of plot you want to make for example...
# below makes a histogram
ggplot(d,aes(bbch.v))+geom_histogram()
### to group object by a certain factor in the data asign an aes() in the geom_object
ggplot(d,aes(bbch.v))+geom_histogram(aes(fill=specificEpithet))
#now change it to geom_histogram(aes(),fill="skyblue")
#what happens when you put something like fill, color or shape outside of the aes() vs inside of it?
?ggplot() # this is how you ask R for help about a given function.

### task 5a: remake the historgram above with each species side by side
# and use the flower bbch values

##### Make a bar plot show how many fruits their are perspecies
# this might help ggplot()+geom_bar()+facet_wrap(~specificEpithet)
 

# Task 6: Now make the heat map from last week (bbch.v vs bbch.f) using
# the geom stat_bin2d() 