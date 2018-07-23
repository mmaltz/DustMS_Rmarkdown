####################################
## dust data
## July 04 2018 
## CN for MM
####################################
## goals - read in mapping file, OTU tables, and distance matrices (weighted or unweighted)  
# assess community composition differences between Dust treatment groups (possibly also by Month species and the interaction), 
# work with denovo OTU table which used numeric values, (as opposed to presence absence) or percentages; 

# expand stats using distance matrices. 
# PERMANOVA,exploring other community analyses, RDA

## preliminaries 
library(dplyr)
## preliminaries 

library(dplyr) ## for data wrangling - %>% function
library(reshape2) ##melt and cast data
library(tidyr) # 'separate' function
library(readxl) #read xlsx files into r on mac computer
library(vegan) # dissimilarity matrix, permanova functions


## data
# Two years
# Four elevations
# Two gene regions: 16S and ITS

# using files in RdataDust
setwd('/Users/maltz/Desktop/RdataDust')#downlaoded location of files
list.files()#show filenames in current working directory
# it is possible to get these directly from google drive into R without downloading...can show you this if that is useful for you

####################################
## browsing computer files from in R, reading in, trouleshooting
# familiarize with file paths & types, searching, pitfalls of data import, indexing, lists, 

list.files()#show filenames in current working directory
list.files() #files in 'Dust' folder in the working directory  
list.files()[3] # third file in that folder (for me it's"")
list.files() # and so on
#Why are there no txt files in the 'RdataDust')
list.files( ,pattern = '.txt') #all txt files

list.dirs()#everything in folder 

# reading in similar file types from different folders
# first - get unique pattern for files
# make a list of all the files
# read them all in

# read in mapping file
list.files()
# MAPPING > SierraMap6.txt
# column H in the mapping file, ‘TreatName’ hmnas of the samples as given on the OTU table in the otu folder,
#library(readxl)
##I got the mapping file to read in well, and named that object 'map'
map<-map
map1<-read_excel("SierraMap6.xlsx", col_names = c('SampleID','BarcodeSequence','LinkerPrimer','Year','Month','SiteCode','RepNum','SiteRep','Site','Elevation','DateCode','DescName','Description'))

# Now I am messing around to try and get the workflow to reflect the best files to read in and the best way to format the columns
map6<-read.csv("SierraMap6.csv", col_names = c('SampleID','BarcodeSequence','LinkerPrimer','Year','Month','SiteCode','RepNum','SiteRep','Site','Elevation','DateCode','DescName','Description'))

map6<-read.csv("SierraMap6.csv", header = TRUE)
View(map6)

#'#SampleID	BarcodeSequence	LinkerPrimerSequence	Year	Month	SiteCode	RepNum	SiteRep	Site	Elevation	DateCode	DescName	Description
#'
View(map6)

# read in otu tables  for each gene group, and the functional group from FUNGuild
# My working directory > RdataDust > filtered_table_w_metadata.txt
list.files()
#head.list<-read.table("headlist.txt", fill=TRUE, stringsAsFactors = FALSE)
# list of the folders that contain report files
filt<-list.files(pattern ='filtered', full.names=TRUE)
#filt
#list.files(filt) #list the files in those folders
# can see there are same files within each folder
#list.files(filt, pattern='table')

# generate list of the file path to each otu table
#This is not working, nor is it necessary -- only one workign directory (only one folder) with one table: otu path -- only needs to read in filtered_table_w_metadata.txt (or .csv or .xlsx)
###But, the function wants to use 'otu_path' and the function wants to use 'folders'; How to modify function to work with the one file that I need to read in, but still modify it effectively?
###having trouble using lapply, and also header_fun and also header_melt_fun
otu_path<-paste0(, filtered_table_w_metadata.txt)
otu_path # should be the otu file paths

## reading in files & common error troubleshooting 
# read in otu files
#library(readxl) #read xlsx files into r on mac computer
# guide: https://github.com/tidyverse/readxl

# if you only have a few files, it is easy enough to read them in individually

# read in single file of '.xls' format
#d16s<-read_excel("filtered_table_w_metadata.xlsx")
A16s<-read.csv("tmp.csv", header = T, stringsAsFactors = F)
View(A16s)
#View(d16s)
#Table had a header with taxonomy as its own column

#du16s<-read.table("filtered_table_w_metadata.txt")
#View(du16s)
# Error:line 1 did not have 17 elements
# usually means that column headers don't match up with the number of columns

# how many columns are in each row? 
count.fields("filtered_table_w_metadata.txt") # number of cols in each row are not equal
#### this is b/c delimination is not clear - in this case b/c cells are tabbed &
#### R is reading in the 'taxonomy' column to have different lengths based on resolution

?read.table
# fix - use 'fill' to fill in missing cells so all rows are equal
du16s<-read.table("filtered_table_w_metadata.txt", fill=TRUE)
View(du16s) # but the first line is the column header
du16s<-read.table("filtered_table_w_metadata.txt", fill=TRUE, header=TRUE)
# when rows are equal, this designates the column names from the first row, but not here

# instead, skip line 1 and assign own names
#du16s<-read.table("filtered_table_w_metadata.txt", fill=TRUE, skip=1, 
#                col.names=c('OJ514','JJ415','JT115','JJ515','OJ615','OJ415','JP414','AP114','JT314','OJ515','JJ114','JT415','JS115','OP515','JP315','OP315','OP615','OS615','OS115','OT515','OT415','AS114','OT414','AT414','OT115','AJ414','OS415','JS415','OS214','JP515','OP614',
#                           'domain','phylum','class','order','family','genus','species','epithet'))
##OTUID not lined up with that column, why? Do I need to type that in?
#View(du16s)# great...but what if it is too tedious to type out all the column names?

dust16s<-read.table("filtered_table_w_metadata.txt", fill=TRUE, skip=1, 
                    col.names=c('OTUID','OJ514','JJ415','JT115','JJ515','OJ615','OJ415','JP414','AP114','JT314','OJ515','JJ114','JT415','JS115','OP515','JP315','OP315','OP615','OS615','OS115','OT515','OT415','AS114','OT414','AT414','OT115','AJ414','OS415','JS415','OS214','JP515','OP614',
                                'domain','phylum','class','order','family','genus','species','epithet'))
View(dust16s)
# great...but what if it is too tedious to type out all the column names?
#####dust16s looks good~!
#
#list.files()#show filenames in current working directory
# read in mapping file
mapSierra<-read_excel("SierraMap6.xlsx")
#View(mapSierra)
#str(mapSierra)
map.meta<-(mapSierra)
View(map.meta)
###
#can save the column names from the file itself
#heads<-read.table("filtered_table_w_metadata.txt", fill=TRUE, stringsAsFactors = FALSE)
#heads<-read.csv("filtered_table_w_metadata.csv", fill=TRUE, stringsAsFactors = FALSE)
# using brackets:
# df = [rows,columns]
# can subset any data frame by indicating the rows and cols wanted 
#heads[1,] # first row containing the column names
##heads[0,]
#heads[1:32,] # rows 1-9
#heads[c(1,2,4,8,9:13),] # rows 1, 2, 4, 8, 9, 10, 11, 12, 13
#heads[1:5,1:5]
#heads[1,1:32]
#heads[1:3,1:32]
#heads[1:5, -9] #all columns except the 9th
#heads[1,1:32] # the given headers are in cols 1:110
#heads[1,-c(33:38)] #or all cols but 12:17  - same
#head.list<-as.character(heads[1,1:32]) #save row and convert to character string
#head.list# using only 1:32 above because we will replace 'taxonomy' with 'domain'



###2018/07/19 head.list isn't working the way I want it to, it is doing the header for the first line entry, not the column headers. Why?
#read.table("headlist.txt")
#View(head.list)
# assign the column names to be the given headers + taxonomy 
#taxa_list<-c('domain','phylum','class','order','family','genus','species','epithet') #taxonomy column headers
#B16<-read.table("filtered_table_w_metadata.txt", fill=TRUE, skip=1, 
#               col.names=c(head.list, taxa_list))
#View(B16)

#B16<-read.table("filtered_table_w_metadata.txt", fill=TRUE,
#               col.names=c(head.list, taxa_list))
#View(B16)
###


# can save the column names from the file itself
#heads<-read.table("filtered_table_w_metadata.txt", fill=TRUE, stringsAsFactors = FALSE)
#heads<-read.csv("filtered_table_w_metadata.csv", fill=TRUE, stringsAsFactors = FALSE)
#heads<-read_excel("filtered_table_w_metadata.xlsx", fill=TRUE, stringsAsFactors = FALSE)
# using brackets:
# df = [rows,columns]
# can subset any data frame by indicating the rows and cols wanted 
#heads[1,] # first row containing the column names
#heads[1:9,] # rows 1-9
#heads[c(1,2,4,8,9:13),] # rows 1, 2, 4, 8, 9, 10, 11, 12, 13
#heads[1:5,1:5]
#heads[1:5, -9] #all columns except the 9th
#heads[1,1:11] # the given headers are in cols 1:11
#heads[1,-c(33:39)] #or all cols but 12:17  - same
#head.list<-as.character(heads[1,1:32]) #save row and convert to character string
#head.list# using only 1:10 above because we will replace 'taxonomy' with 'kingdom'

# assign the column names to be the given headers + taxonomy 
#taxa_list<-c('domain','phylum','class','order','family','genus','species') #taxonomy column headers
#Dust16s<-read.table("filtered_table_w_metadata.txt", fill=TRUE, skip=1, 
#                col.names=c(head.list, taxa_list))
##still got error: more columns than column names
View(dust16s)
# especially helpful to automate work pipelines and combine multiple data sheets w/ different column names

####################################
## functions, reshaping and combining dataframes, 'apply', automating file imports
# writing a function to read in multiple files, rename columns, and combine 
# when they are same type of data - i.e. otu tables for each Month & gene region are separate files

## goal - read in all 4 otu tables and combine into single df for multivariate tests
#otu_path
#otu_path[2] # select the 2nd element - when it's a list, here is only 1 dimension in brackets 

# a function to read in otu files and add taxonomy as last 7 columns
################@@@@@@@@
#header_fun<-function(x){
#  taxa_list<-c('kingdom','phylum','class','order','family','genus','species') #taxonomy column headers
#  temp.df<-read.table(otu_path[1], fill=TRUE, stringsAsFactors = FALSE)
#  head_list<-as.character(temp.df[1, 1:(length(temp.df)-7)]) #dataframe headers
#  out.df<-read.table(x, fill=TRUE, skip=1,col.names=c(head_list,taxa_list)) # read in and assign headers
#  out.df$file<-paste(x) # create additional column that designates source file
#  return(out.df)
#}

#@#######Experiment
#header_fun<-function(x)
#  {
#  taxa_list<-c('kingdom','phylum','class','order','family','genus','species') #taxonomy column headers
#View(taxa_list)  
#  temp.df<-read.table(otu_path[1], fill=TRUE, stringsAsFactors = FALSE)
#View(temp.df)  
#  head_list<-as.character(temp.df[1, 1:(length(temp.df)-7)]) #dataframe headers
#View(head_list)  
#  out.df<-read.table("RdataDust/Dust16s_report/otu/otu_taxon.xls", fill=TRUE, skip=1,col.names=c(head_list,taxa_list)) # read in and assign headers
#View(out.df)  
#  out.df$file<-paste("RdataDust/Dust16s_report/otu/otu_taxon.xls") # create additional column that designates source file
# View(out.df) 
#  return(out.df)
#}  

# file1<-
#}
?return

###################
header_fun<-function(x){
  taxa_list<-c('kingdom','phylum','class','order','family','genus','species') #taxonomy column headers
  temp.df<-read.table(otu_path[1], fill=TRUE, stringsAsFactors = FALSE)
  head_list<-as.character(temp.df[1, 1:(length(temp.df)-7)]) #dataframe headers
  out.df<-read.table(x, fill=TRUE, skip=1,col.names=c(head_list,taxa_list)) # read in and assign headers
  out.df$file<-paste(x) # create additional column that designates source file
  return(out.df)
}  
###################


#######################################
#Deconstructing the function_header_fun
#######################################
taxa_list<-c('domain','phylum','class','order','family','genus','species','epithet')
View(taxa_list)
temp.df<-read.csv("all_otusBD.csv", fill=TRUE, stringsAsFactors = FALSE)
View(temp.df)
head_list<-as.character(temp.df[1, 1:(length(temp.df)-7)]) #dataframe headers
head_list<-as.character(temp.df[1, 1:33]) #dataframe headers

###########$$$$$$$$$$$$$$$$##########
#Issue with head_list)
###########$$$$$$$$$$$$$$$$##########

View(head_list)
#View(head.list)
#now try this function for a single file path (x)
#file1<-header_fun(otu_path[1])
#View(file1)
# other files
#file2<-header_fun(otu_path[2])
#View(file2)
#file3<-header_fun(otu_path[3])
#View(file3)
#file4<-header_fun(otu_path[4])
#View(file4)
#all_otus<-rbind(file1, file2, file3, file4)
#View(all_otus)# now they are all in a single file

# read in all files at once
# 'apply' to a 'list' = 'lapply'
# apply 'modified 'read.table''header_fun' function to the list of otu files (otu_path)
# guide: https://nicercode.github.io/guides/repeating-things/
###########$$$$$$$$$$$$$$$$##########
#Issue with lapply (why isn't this working? Is is because otu_path is moot? Or because header_fun won't work?))
###########$$$$$$$$$$$$$$$$##########
ldf<-lapply(otu_path, header_fun)
str(ldf)# this produces a list of dataframes
length(ldf) # length is 4 because there were 4 files imported, and a dataframe made for each file
#View(ldf[[2]]) # view second df in list

#library(dplyr)
#all_otus<-bind_rows(ldf)##combine list of dataframes
View(dust16s)
all_otus<-dust16s    
View(all_otus)# now they are all in a single file - same result as before

#save
write.csv(all_otus, 'all_otus_16SDust_20180719.csv', row.names=FALSE)
write.csv(all_otus, 'all_otusBD2.csv', row.names=FALSE)

# now this last part (binding the dataframes) will only work if the column headers are exactly the same in each file
# so if there were different numbers of samples between the files...

# melt dataframe so there is a column for each sample name
# basically reduce a subset of columns into 2 columns (1 with the former column name, 1 with the values)
#library(reshape2)
View(all_otus)

###########$$$$$$$$$$$$$$$$##########
#Issue with file2)
###########$$$$$$$$$$$$$$$$##########
file2.melt<-melt(all_otus, id.vars=c('OTUID',taxa_list))
View(file2.melt) #compare this to file2

# include syntax to name new columns
file2.melt<-melt(all_otus, id.vars=c('OTUID',taxa_list), value.name='count', variable.name='sample')
str(file2.melt)
View(file2.melt)
all_otus_melt<-file2.melt
View(all_otus_melt)
# to consider another time - 
# how to deal with otus that are same taxomony but different otu?
# are the otu id's 'denovo0' the same otu across files? or these are unique to each file?

## add 'melt' to function
header_melt_fun<-function(all_otus){
  taxa_list<-c('domain','phylum','class','order','family','genus','species') #taxonomy column headers
  temp.df<-read.table(all_otus, fill=TRUE, stringsAsFactors = FALSE)
  head_list<-as.character(temp.df[1, 1:(length(temp.df)-7)]) #dataframe headers
  out.df<-read.table(all_otus, fill=TRUE, skip=1,col.names=c(head_list,taxa_list)) # read in and assign headers
  out.df$file<-paste(all_otus) # create additional column that designates source file
  out.df.melt<-melt(out.df, id.vars=c('ID','file', taxa_list), value.name='count', variable.name='sample')
  return(out.df.melt)
}
View(header_melt_fun)
#####
#####This isn't working'
###########$$$$$$$$$$$$$$$$##########
#Issue with head_melt_fun
###########$$$$$$$$$$$$$$$$##########
###########$$$$$$$$$$$$$$$$##########
#Issue with lapply)
###########$$$$$$$$$$$$$$$$##########
all_otus_melt<-bind_rows(lapply(all_otus, header_melt_fun))
View(all_otus_melt) #everything's here

write.csv(all_otus_melt, 'all_otus_melt20180719.csv', row.names=FALSE)



####################################
## Dust data
## July19 2018 
##ModifydustDust
####################################
####################################
## dust data
## July 18 2018 
## CN for MM
####################################
## goals - read in mapping file, OTU tables, and distance matrices (weighted or unweighted)  
# assess community composition differences between Dust treatment groups (possibly also by Month species and the interaction), 
# work with denovo OTU table which used numeric values, (as opposed to presence absence) or percentages; 

# expand stats using distance matrices. 
# PERMANOVA,exploring other community analyses, RDA

## preliminaries 
library(dplyr)
## preliminaries 

library(dplyr) ## for data wrangling - %>% function
library(reshape2) ##melt and cast data
library(tidyr) # 'separate' function
library(readxl) #read xlsx files into r on mac computer
library(vegan) # dissimilarity matrix, permanova functions


## data
# Two years
# Four elevations
# Two gene regions: 16S and ITS

# using files in RdataDust
setwd('/Users/maltz/Desktop/RdataDust')#downlaoded location of files
list.files()#show filenames in current working directory

####################################
#Modify dust SCRIPTS
####################################
####################################
## for Dust 16S
####################################
###################################
####################################
#R CLASS
######################################
## Part 2 - data wrangling 
# By CN
# 04/03/2018
#'/Users/maltz/Desktop/RdataDust/
# Preliminaries
setwd('/Users/maltz/Desktop/RdataDust')

list.files()
library(readr)
all_otus_melt <- read_csv("all_otus_melt.csv")
View(all_otus_melt)
library(readxl)
library(dplyr)
library(reshape2)
library(tidyr)
library(vegan)
library(ggplot2)# ggplot resource -http://rpubs.com/collnell/ggplot2

library(tidyverse) # useful packages in 1 - dplyr, ggplot2, tidyr +


######################################
## creating statistical pipelines 
# making workflow and analyses reproducbile, efficient, easy to repurpose

## building pipelines - think of everything in pieces, checking at each step to make sure has desired effect
# then, reconstruct working pieces into short & understandable chunks, clean it up
# often not obvious where steps can be simplified or piped together during first work through, ordering
# if code is not working - put '#' after lines to locate problem 

## remember - '%>%" is called a pipe or a chain
# it chains together a series of commands administered to the same data frame
# common across 'tidyverse' packages and many more
# benefit - use consistent syntax, organized and easy to read, do many things in one chunk, reduces typing & related errors

## cleaning up progress made last week (dust_wrangling.R) 
# streamline data cleaning workflow

## read in data and clean up variable codings
library(data.table) #fread function
list.files()

#otus<-fread('all_otus.melt.csv')%>% # reads in large datasets faster than 'read.csv'
# filter(count >= 1)%>% # drop 0 abundance data
#unite(col='taxa', kingdom:species, remove=TRUE)%>%
#separate(col='file', into=c('Elevation','gene'), sep='_')%>% # separate into separate columns for Elevation and gene
#group_by(Month, gene, sample, taxa, ID)%>% # set grouping level you want data aggregated at
#summarize(otu_abun = sum(count), otu_rich =length(unique(ID))) #%>% # summarize the abundances at desired level - here is just each unique ID x taxa (so every otu)

otus<-fread('all_otus_melt.csv')%>% # reads in large datasets faster than 'read.csv'
  filter(count >= 1)%>% # drop 0 abundance data
  unite(col='taxa', domain:epithet, remove=TRUE)%>%
  group_by(sample, taxa)%>% # set grouping level you want data aggregated at
  summarize(otu_abun = sum(count), otu_rich =length(unique(OTUID))) #%>% # summarize the abundances at desired level - here is just each unique ID x taxa (so every otu)

## otu abundance vs otu richness
# otu_abun =e total abundance (counting multiple reads of the same otu), otu_rich = # of unique IDs - denovos
# here, otu_rich should be 1 because we are summarizing at the ID level, BUT if done differently, will tell you the richness of otus vs the abundance
# worth considerating - what is more meaningful, total reads or unique otus? are reads skewed by methodology or do they reflect abundance in an accurate way?
# also - how to deal with multiple 'ID's of same 'taxa'? read as a distinct otu, but taxanomic difference unknown - consistent across Month species

head(otus) # in long format, where you can easily bind with alternative grouping data (functional groups)
write.csv(otus, '/Users/maltz/Desktop/RdataDust/otus_by_OTUID20180719.csv', row.names=FALSE)

## make community data matrix for analysis
# summarize data at desired groupings for community data matrix (what do you want your 'species' or columns to be?)
## If you are going to regroup that into another taxonomic split -- species  
# summarize at the class level to look at general patterns in community composition
# for other - can use 'left_join' to bind data, then do these next steps

otu.class<-otus%>%
  separate(taxa, into=c('domain','phylum','class','order','family','genus','species','epithet'), sep=';')%>% #this produces a warning 'too few values' because not all taxa are id'd to species, not an issue
  dplyr::select(class, sample, otu_abun, otu_rich)%>% # drop other taxanomic columns not using
  mutate(class = tolower(as.character(class)))%>%
  mutate(class = gsub('_c', '', class))%>% # clean up class names for consistency
  mutate(class = gsub('_', '', class))%>% 
  group_by(sample, class)%>% # the lowest grouping level here (class) is what you want your columns to be
  summarize(otu_abun= sum(otu_abun), otu_rich = sum(otu_rich))#%>% # summarize at the class level
head(otu.class) # df that tells you the otu abundance and richness for each sample of each chamber status at the class level (long format)
write.csv(otu.class, '/Users/maltz/Desktop/RdataDust/otu_classBD.csv', row.names=FALSE)
# note - different than before I am using the 'mutate' function with 'tolower', 'gsub' and 'recode' within them
# this works b/c mutate makes new variables, and those other commands work on vectors, not dataframes
# makes it easier to pipe together for flow and don't need to write 'otu.class$' in front of everything

unique(otu.class$class)# use to verify class names - duplicates? misspellings? NA are those that didn't get id'd to class
###################$$$$$$$$$$$$$$$$$$$$$$$$$$$$###############
#Change the names of the 'hloroplast' to 'chloroplast'
###################$$$$$$$$$$$$$$#############################
## cast into community dataframe (samples as rows, species as columns), & add grouping variables (WIDE DATA)
# here need to decide whether total otu abundance or the richness of reads matters more - use to cast
otu.cast<-otu.class%>%
  dcast(sample~class, value.var='otu_abun')#%>% # for each sample (lowest level within sample), otu abundance of each class (as columns)
# replace NA with 0 in abundances
otu.cast[is.na(otu.cast)]<- 0 
head(otu.cast)

## read in mapping data
# remember - not all sample and treatnames matched due to an R in some names
#map.meta<-read_excel("mapBr2.xlsx")
# dplyr::select(-contains('Sequence'), -Variable, -Description, -`#SampleID`)%>% # drop cols not needed
# mutate(treat = recode(TreatName, "E15010"='E150R10', "E15011"="E150R11", "E15012"="E150R12", .default = TreatName)) # default gives variable value if not specifically named in 'recode'
View(map.meta)
unique(map.meta$DescName)
colnames(map.meta)
write.csv(map.meta, '/Users/maltz/Desktop/RdataDust/BD_metadata20180719.csv', row.names=FALSE)# rewite so don't need to recode variables again
map_meta<-read.csv("BD_metadata1.csv")
# combine
otu_map<-left_join(otu.cast, map_meta, by=c('sample'))%>%
  dplyr::select(sample,Year,Month,SiteCode,RepNum,SiteRep,Site,Elevation,DateCode,Description,everything())
View(otu_map)

##This only has class
write.csv(otu_map, 'otus_by_classBD.csv', row.names=FALSE)
#View


##############################################################################################################
######################################
## Part 3 - multivariate community analyses and data visualization
# By CN
# 02/07/2018

# Preliminaries
#dust_otu.txt
#mapBr2.xlsx
setwd('/Users/maltz/Desktop/RdataDust')
library(dplyr)
library(reshape2)
library(tidyr)
library(vegan)
library(ggplot2)# ggplot resource -http://rpubs.com/collnell/ggplot2

library(tidyverse) # useful packages in 1 - dplyr, ggplot2, tidyr +

# community matrix at order level with mapping data
#otu_map<-read.csv('/Users/maltz/Desktop/RdataDust/otus_by_classBD.csv')
head(otu_map)

colnames(otu_map)
comm.grps<-otu_map%>%dplyr::select(sample:Description) #mapping data
colnames(comm.grps)

comm.mat<-otu_map%>%dplyr::select(-c(sample:Description)) # community matrix - all but mapping data

######################################
## comparing ecological communities
# diversity vs composition
# abundance and richness are univariate response variables used to quantify communities
# in multivariate analyses we have these variables for multiple entities
# similarly, multivariate analyses have counterparts in univariate stats - t-test, ANOVA, mutliple regression

## univariate analyses of diversity

head(comm.grps)
str(comm.grps)
View(comm.grps)
comm.mat$hloroplast<-as.numeric(comm.mat$hloroplast)
comm.mat$ktedonobacteria<-as.numeric(comm.mat$ktedonobacteria)
comm.mat$sphingobacteriia<-as.numeric(comm.mat$sphingobacteriia)
comm.mat$phycisphaerae<-as.numeric(comm.mat$phycisphaerae)
comm.mat$betaproteobacteria<-as.numeric(comm.mat$betaproteobacteria)
comm.mat$bacilli<-as.numeric(comm.mat$bacilli)
comm.mat$deltaproteobacteria<-as.numeric(comm.mat$deltaproteobacteria)
comm.mat$alphaproteobacteria<-as.numeric(comm.mat$alphaproteobacteria)
comm.mat$'[spartobacteria]'<-as.numeric(comm.mat$'[spartobacteria]')
comm.mat$acidimicrobiia<-as.numeric(comm.mat$acidimicrobiia)
comm.mat$'acidobacteria-6' <-as.numeric(comm.mat$'acidobacteria-6')
comm.mat$acidobacteriia<-as.numeric(comm.mat$acidobacteriia)
comm.mat$actinobacteria<-as.numeric(comm.mat$actinobacteria)


comm.mat$acidobacteriia<-as.numeric(comm.mat$acidobacteriia)

comm.mat$thermoleophilia<-as.numeric(comm.mat$thermoleophilia)
comm.mat$thermomicrobia<-as.numeric(comm.mat$thermomicrobia)
comm.mat$thermotogae<-as.numeric(comm.mat$thermotogae)
comm.mat$tk10<-as.numeric(comm.mat$tk10)
comm.mat$tk17<-as.numeric(comm.mat$tk17)
comm.mat$tm1<-as.numeric(comm.mat$tm1)
comm.mat$'tm7-1'<-as.numeric(comm.mat$'tm7-1')

comm.mat$'tm7-3'<-as.numeric(comm.mat$'tm7-3')
comm.mat$vadinha49<-as.numeric(comm.mat$vadinha49)
comm.mat$vc21bac22<-as.numeric(comm.mat$vc21bac22)
comm.mat$'verruco-5'<-as.numeric(comm.mat$'verruco-5')

comm.mat$verrucomicrobiae<-as.numeric(comm.mat$verrucomicrobiae)
comm.mat$planctomycetia<-as.numeric(comm.mat$planctomycetia)
comm.mat$zb2<-as.numeric(comm.mat$zb2)
comm.mat$'NA'<-as.numeric(comm.mat$'NA')
#comm.mat$tk17<-as.numeric(comm.mat$tk17)

str(comm.mat2)
tail(comm.mat2)
View(comm.mat2)
## does diversity vary across groups?
# compute diversity indices
indices <- comm.grps
#indices <- comm.mat
##Error that >0 not working on factors
#comm.mat2$solibacteres<-as.numeric(comm.mat$solibacteres)
indices$richness <- rowSums(comm.mat2>0)
indices$shannon <- diversity(comm.mat2, index='shannon')
indices$rarified <- c(rarefy(comm.mat2, sample=1192)) # rarefied diversity for a given sample size
##Rarefy to a relevent number for this dataset!!!

## visualize differences in diversity by Month species
ggplot(indices, aes(x = Elevation, y = richness))+geom_boxplot() # seemingly higher diversity for Month T

# what about when rarified?
# AND color by Elevation 
ggplot(indices, aes(x = Elevation, y = rarified))+geom_boxplot(aes(fill = elevation))
########################################$$$$$$$$$$$$$$$$$$
#Rarified not found #######################################
########################################$$$$$$$$$$$$$$$$$$$
ggplot(indices, aes(group=SiteCode, x = DateCode, y = richness))+geom_boxplot() # seemingly higher diversity for tree T


# common points of confusion
ggplot(indices, aes(group=SiteCode, x = Elevation, y = rarified))+geom_boxplot(aes(color = Elevation)) # fill vs color - what is modified depends on 'geom' type

ggplot(indices, aes(group=SiteCode, x = Elevation, y = rarified))+geom_boxplot(aes(color = Elevation)) # mapping color in 'aes' vs outside
ggplot(indices, aes(group=SiteCode, x = Elevation, y = rarified))+geom_boxplot(aes(color = 'SiteCode'))
ggplot(indices, aes(group=Month, x = Elevation, y = rarified))+geom_boxplot(color = 'blue') # set specific color to ALL
ggplot(indices, aes(group=SiteCode, x = Elevation, y = rarified))+geom_boxplot(aes(color = Elevation, fill = Elevation)) # within 'aes' maps to the levels of a variable
ggplot(indices, aes(group=SiteCode, x = Elevation, y = rarified))+geom_boxplot(color = 'darkslateblue', fill = 'yellow') 

ggplot(indices, aes(group=SiteCode, x = Elevation, y = rarified))+geom_boxplot(aes(color = Elevation)) # fill vs color - what is modified depends on 'geom' type

ggplot(indices, aes(group=Elevation, x = Elevation, y = rarified))+geom_boxplot(aes(color = Month)) # mapping color in 'aes' vs outside
ggplot(indices, aes(group=SiteCode, x = Elevation, y = rarified))+geom_boxplot(aes(color = 'SiteCode'))
ggplot(indices, aes(group=Month, x = Elevation, y = rarified))+geom_boxplot(color = 'blue') # set specific color to ALL
ggplot(indices, aes(x = Elevation, y = rarified))+geom_boxplot(aes(color = Elevation, fill = Elevation)) # within 'aes' maps to the levels of a variable
ggplot(indices, aes(x = Elevation, y = rarified))+geom_boxplot(color = 'darkslateblue', fill = 'yellow') 


#plot points on top of boxplot
ggplot(indices, aes(group=Elevation, x = Elevation, y = rarified))+
  geom_boxplot(aes(color = Elevation))+
  geom_point(aes(color=Elevation), size=3)

# ANOVA - differences by groups
Elevation.rich.aov<-aov(rarified~Elevation, data=indices)
summary(Elevation.rich.aov)
# order richness is significantly higher for Elevation (low) compared to elevation (high)

# what variables explain richness?
pairs(indices%>%dplyr::select(Elevation, Month, Site, DateCode, SiteRep, richness)) # data exploration
# looks like Elevation species affects TotalP, OrganicMatter, pH
# same with dust

# linear regression - trends
rich.lm<-lm(richness~SiteCode*Elevation, data=indices)
summary(rich.lm)
summary(aov(rich.lm))

#linear regression - trends
richTime.lm<-lm(richness~Month*Elevation, data=indices)
summary(richTime.lm)
summary(aov(richTime.lm))

richTime2.lm<-lm(richness~DateCode*Elevation, data=indices)
summary(richTime2.lm)
summary(aov(richTime2.lm))

richTime3.lm<-lm(richness~Year*Elevation, data=indices)
summary(richTime3.lm)
summary(aov(richTime3.lm))

# for interactions, need to consider SS type
install.packages(car)
??MASS

??ANOVA

library(MASS)

car::Anova(rich.lm)
car::Anova(rich.lm, type='III')

# examine distribution of residuals
resids<-resid(rich.lm)
shapiro.test(resids)
plot(resids)
qqnorm(resids)
qqline(resids)

# include environmental variables
rich.env.lm<-lm(richness~Elevation+Month+DateCode, data=indices)
summary(aov(rich.env.lm))

#mean and sd richness by species
Elevation.rich<-indices%>%group_by(Elevation)%>%summarize(mean = mean(richness), sd=sd(richness))
#summary(aov(Elevation.rich))
?indices
??sd
#plot means
ggplot(Elevation.rich, aes(Elevation, mean))+
  geom_point(size=3)

##error bars
ggplot(Elevation.rich, aes(Elevation, mean))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd))

ggplot(Elevation.rich, aes(Elevation, mean))+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0)+
  geom_point(size=3)

# editing plot themes
ggplot(Elevation.rich, aes(Elevation, mean))+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0)+
  geom_point(size=3)+
  theme(panel.background = element_rect(fill='white'))## change background color

# black axis lines
myplot<-ggplot(Elevation.rich, aes(Elevation, mean))+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0)+
  geom_point(size=3)+
  theme(panel.background = element_rect(fill='white'), axis.line = element_line(color='black'))
myplot

myplot+theme(axis.text = element_text(size=10)) #tick label text size
myplot+theme(axis.text = element_text(size=10), axis.title=element_text(size=12)) #axis label size

# set axis limits
myplot+ylim(0,55)
myplot+ylim(0,NA)

mynewplot<-myplot+geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1)+ylim(0,NA)+
  theme(axis.text = element_text(size=12), axis.title=element_text(size=14))+
  labs(x='Elevation species', y='Richness')
mynewplot

## changing colors
myplot<-ggplot(Elevation.rich, aes(Elevation, mean))+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.1)+
  geom_point(size=3, aes(color=Elevation))+
  theme(panel.background = element_rect(fill='white'), axis.line = element_line(color='black'),
        axis.text = element_text(size=12), axis.title=element_text(size=14))+
  ylim(0,NA)+
  labs(x='Elevation species', y='Richness')
myplot

myplot+theme(legend.position='bottom')
myplot+theme(legend.position='left')
myplot+theme(legend.position='none')

myplot+theme(legend.position='none')+
  scale_color_manual(values=c('purple','green'))
# can use anyvalid names - https://www.nceas.ucsb.edu/~frazier/RSpatialGuides/colorPaletteCheatsheet.pdf

myplot+theme(legend.position='none')+
  scale_color_manual(values=c('purple','green'))
# any scale you set (shape, color, fill, etc) you can change with a similar syntax
# i.e. scale_fill_manual(), scale_alpha_manual

######################################
## multivariate analyses
# similar tests depending on type of data we have, how much, nature of distribution
# ordination, dimension reduction, gradient analyses

## general questions
# Do groups differ in composition? -manova/permanova, NMDS
# How do the data group naturally, what is most similar? - cluster analysis, 
# What best explains/predicts variation in composition between groups? - RDA, dbRDA, discriminant analysis, simper, random forest/CART
# How do multiple explanatory variables relate to multiple outcomes? - canoncial correlation analysis

## data exploration - correlations between variables
install.packages('corrplot')
library(corrplot)

##correlations - can only use numeric variables
comm.num<-otu_map%>% dplyr::select_if(is.numeric)
str(comm.num)

# look at correlations
comm.cor<-cor(comm.num)
corrplot(comm.cor) # hard to digest
corrplot(comm.cor, order='alphabet', tl.col='black')

# just classes
corrplot(cor(comm.mat), order='hclust', hclust.method='average', tl.col='black', tl.cex=.75)
# some structure in data - blocks indicte co-occuring groups
# i.e. entrrhizales:hltermanniales co-occuring & micrascales:atheliales

## first - do groups differ? permanova
comm.rel<-decostand(comm.mat2, method='total')# relative composition
comm.dis<-vegdist(comm.mat2, method='bray')# dissimialrity

## permanova - multivariate ANOVA
Elevation.perm<-adonis(comm.dis~Elevation, data=comm.grps, permutations=10000)
Elevation.perm

DateCode.Elevation.perm<-adonis(comm.dis~Elevation*DateCode, data=comm.grps, permutations=10000)
DateCode.Elevation.perm

Month.Elevation.perm<-adonis(comm.dis~Elevation*Month, data=comm.grps, permutations=10000)
Month.Elevation.perm

## include random effects or blocking
# using 'strata' in 'adonis2' restricts the permutations within that variable
# like a blocking variable or random effect for the permanova
dust.perm<-adonis2(comm.dis~Month, data=comm.grps, permutations=10000, strata='Elevation')
dust.perm # so even though overall diversity did not differ with dust level, there are compositional changes

## investigate similarity among samples

## hierarchical clustering
comm.hclust<-hclust(comm.dis, method='average')
plot(comm.hclust, labels=comm.grps$sample) # numbers correspond to row numbers unless declared
# grouped based on compositional similarity
plot(comm.hclust, labels=comm.grps$Month)
plot(comm.hclust, labels=comm.grps$SiteRep)
plot(comm.hclust, labels=comm.grps$SiteCode)
plot(comm.hclust, labels=comm.grps$Year)


## NMDS ordination 
# nmds - ideal for community data due to prevalence of rare species, non-euclidean data (count data)
comm.nmds<-metaMDS(comm.dis)
comm.nmds<-metaMDS(comm.rel, distance='bray') # can give either community matrix of ready made dissimilarities

# depdning on data transformationa nd distance method, mds config and stress may be poor
metaMDS(comm.mat, distance='bray') # applies sqrt and wisconsin double standardization, high stress
# likely not a good visual representation of data- try other

# stress
comm.nmds #  stress = 0.13 good
stressplot(comm.nmds) # observed vs cnonfigured distances. closer to R2 =1 is ideal
# basically shows - does the mds configuration represent the true dissimilarities among our data smamples?
# product of rank-ordered distance method of nmds
ordiplot(comm.nmds, display='sites', type='text')
ordipointlabel(comm.nmds) # hard to maipulate this to look as you want

# overlay the cluster diagram above
ordiellipse(comm.nmds, comm.grps$Elevation, conf=0.95, label=TRUE) # ovelray 95%CI
ordicluster(comm.nmds, comm.hclust, col='gray')# clsuter - hard to see

## to plot nmds is ggplot2:
# extract NMDS site scores (samples)
nmds.sites<-as.data.frame(comm.nmds$points)
View(nmds.sites)# note - rows correspond to the samples 

# there are a variaty of ways you can recreate this dataframe if row ordering is consistent
#manually
nmds.sites$Elevation<-comm.grps$Elevation
nmds.sites$Month<-comm.grps$Month
nmds.sites$Year<-comm.grps$Year
nmds.sites$SiteCode<-comm.grps$SiteCode
# and so on...

#####################################
#Error in this command
#####################################

#or can implement in the creation of anew df
comm.nmds$points[,2]#the first column MDS1

nmds.sites<-data.frame(NMDS1 = comm.nmds$points[,1],
                       NMDS2 = comm.nmds$points[,2],
                       Elevation = comm.grps$Elevation,
                       rep = comm.grps$level,
                       group = comm.grps$pH,
                       SID = comm.grps$OrganicMatter)

## plot samples in NMDS
ggplot(data=nmds.sites)+
  geom_point(aes(x=NMDS1, y=NMDS2))

## explore how composition is grouping based on mapping variables
#color by Elevation species
ggplot(data=nmds.sites)+
  geom_point(aes(x=NMDS1, y=NMDS2, shape=Elevation), size=3)
#change shapes for dust level
ggplot(data=nmds.sites)+
  geom_point(aes(x=NMDS1, y=NMDS2, shape=Elevation, color=level), size=3)
#cahgne size by pH
ggplot(data=nmds.sites)+
  geom_point(aes(x=NMDS1, y=NMDS2, shape=Month, color=Year, size=pH))
# samples generally grouping by dust level, pH?

## species scores - tell us the associations between orders/groups
# extract NMDS species scores for plotting
nmds.sp<-data.frame(NMDS1 = comm.nmds$species[,1],
                    NMDS2 = comm.nmds$species[,2])
View(nmds.sp)
nmds.sp$order<-rownames(nmds.sp) #make variable for names

# add to plot
ggplot(data=nmds.sites)+
  geom_point(aes(x=NMDS1, y=NMDS2, shape=Month, color=level), size=3)+
  geom_text(data=nmds.sp, aes(x=NMDS1, y=NMDS2, label=order))
# lots of technical errors can happen when plotting from different dataframes in same plot
# need to be careful where 'data' is declared, and what is in the 'aes'
# when only one data source, 'data' can be given within the 'ggplot()' to be applied to all other 'geoms'
# unless another data is given in a geom - will apply to all geoms unless otherwise given

# if multiple sources -  make sure data is given for each geom layer
# here I move the data from 'ggplot' to 'geom_point'
# equivalent of plot above 
ggplot()+
  geom_point(data=nmds.sites, aes(x=NMDS1, y=NMDS2, shape=Elevation, color=level), size=3)+
  geom_text(data=nmds.sp, aes(x=NMDS1, y=NMDS2, label=order))
str(nmds.sp)

## there is a certain art to building ordinations - work within constraints of your data to best reflect relationships
## this is a ton of orders - let's filter to only orders than made up > 1% of all otu reads
colSums(comm.mat2) # total in each column
sum(colSums(comm.mat2)) ##total 
100*colSums(comm.mat2)/sum(colSums(comm.mat2)) # percentage of total
100*colSums(comm.mat2)/sum(colSums(comm.mat2)) > 1 # which columns are greater than .1%?

top.classes<-comm.mat2[,100*colSums(comm.mat2)/sum(colSums(comm.mat2)) > .1] # selecct only the columns that are >1
dim(top.classes) ##reduces to 25 orders
# there are a lot of other ways to select variables depending on the hypothesis 
# some filter like this prior to calculating distances if too many rare groups
top.list<-colnames(top.classes)
nmds.sp.top<-nmds.sp%>%filter(order %in% top.list)

## replot a subset of the orders and relationship to samples
ggplot()+
  geom_point(data=nmds.sites, aes(x=NMDS1, y=NMDS2, shape=Month, color=level), size=3)+
  geom_text(data=nmds.sp.top, aes(x=NMDS1, y=NMDS2, label=order))
# exploratory 


# adding to the same plot form before
nmds.plot<-ggplot(data=nmds.sites, aes(x=NMDS1, y=NMDS2))+
  geom_point(size=3, aes(shape=Elevation, color=level))+ 
  geom_text(data=nmds.sp.top, aes(x=NMDS1, y=NMDS2, label=order))+
  theme_minimal()
nmds.plot

nmds.plot+stat_ellipse(aes(color=Elevation, lty=Month)) # lty = linetype


################################################################################################
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ################################################################################################  
## goals - read in mapping file, OTU tables, and distance matrices (weighted or unweighted)  
# assess community composition differences between chamber treatment groups
# work with OTU table which used numeric values, (as opposed to presence absence) or percentages; 
# expand stats using distance matrices. 
# PERMANOVA,exploring other community analyses, RDA
## preliminaries 
library(dplyr)
## data
# Two years
##Four elevations
# using files in RdataDust folder (contianing dust)
setwd('/Users/maltz/Desktop/RdataDust')#downlaoded location of files
list.files()#show filenames in current working directory
list.dirs()#everything in folder 

# read in mapping file
list.files()
# column A in the mapping file, ‘sample’ hmnas of the samples as given on the OTU table,
#map<-read.table("SierraMap6.txt", col.names =c('sample','Elevation','group','Rep','d'))
#map<-read.table("mapSP9.txt")
View(map)

# read in otu table
list.files()
## reading in file & common error troubleshooting 
library(readxl) #read xlsx files into r on mac computer
# guide: https://github.com/tidyverse/readxl
# if you only have a few files, it is easy enough to read them in individually
# read in single file of '.xls' format
dust_otu<-read_excel("dust_otu.xlsx")
count.fields(dust_otu)
#Dotu<-read_excel("dust_otu.xlsx")
# not working for some reason - try as table
#b_otu<-read.table("dust_otu.txt")
# Error:line 43 did not have 22 elements
# usually means that column headers don't match up with the number of columns

# how many columns are in each row? 
count.fields("dust_otu.txt") # number of cols in each row are not equal
# this is b/c delimination is not clear - in this case b/c cells are tabbed &
# R is reading in the 'taxonomy' column to have different lengths based on resolution

?read.table
# fix - use 'fill' to fill in missing cells so all rows are equal
#b_otu<-read.table("dust_otu.txt", fill=TRUE)
#View(b_otu) # but the first line is the column header
dust_otu.txt<-read.table("dust_otu.txt", fill=TRUE, header=TRUE)
# when rows are equal, this designates the column names from the first row, but not here
View(dust_otu.txt)
# instead, skip line 1 and assign own names
dust_otu<-read.table("dust_otu.txt", fill=TRUE, skip=1, 
                        col.names=c('SID','MMCH5','MMCH5R','MMCH1R','MMCH1','MMN4','MMN3R','MMN4R','MMCH3R','MMCH4','MMN3','MMN2','MMCH3','MMN2R','MMCH4R',
                                    'kingdom','phylum','class','order','family','genus','species'))
View(dust_otu)# great...but what if it is too tedious to type out all the column names?

View(taxa_list)

all_otus<-dust_otu
View(all_otus)# now they are all in a single file

library(dplyr)
write.csv(all_otus, '/Users/maltz/Desktop/RdataDust/all_otus4.csv', row.names=FALSE)
list.files()
# now this last part (binding the dataframes) will only work if the column headers are exactly the same in each file
# so if there were different numbers of samples between the files...

# melt dataframe so there is a column for each sample name
# basically reduce a subset of columns into 2 columns (1 with the former column name, 1 with the values)
library(reshape2)
#View(file2)
View(all_otus)
all_otus.melt<-melt(all_otus, id.vars=c('SID',taxa_list))
View(all_otus.melt) 

# include syntax to name new columns
all_otus.melt<-melt(all_otus, id.vars=c('SID',taxa_list), value.name='count', variable.name='sample')
str(all_otus.melt)
# to consider another time - 
# how to deal with otus that are same taxomony but different otu?
# are the otu id's 'denovo0' the same otu across files? or these are unique to each file?

## add 'melt' to function
#header_melt_fun<-function(x){
# taxa_list<-c('kingdom','phylum','class','order','family','genus','species') #taxonomy column headers
#temp.df<-read.table(otu_path[1], fill=TRUE, stringsAsFactors = FALSE)
#head_list<-as.character(temp.df[1, 1:(length(temp.df)-7)]) #dataframe headers
#out.df<-read.table(x, fill=TRUE, skip=1,col.names=c(head_list,taxa_list)) # read in and assign headers
#out.df$file<-paste(x) # create additional column that designates source file
#out.df.melt<-melt(out.df, id.vars=c('ID','file', taxa_list), value.name='count', variable.name='sample')
#return(out.df.melt)
#}

#all_otus_melt2<-bind_rows(lapply(otu_path, header_melt_fun))
View(all_otus.melt) #everything's here

write.csv(all_otus_melt, '/Users/maltz/Desktop/RdataDust/all_otus_melt4.csv', row.names=FALSE)

######################################
## wranging skills
# manipulating data with dplyr, tidyr, and reshape2, working with characeter strings
# familiarity with different functions, splitting and combinging existing variables

setwd('/Users/maltz/Desktop/RdataDust')# downlaoded location of files
list.files()

library(dplyr)
library(reshape2)

# from otus created above, we can regroup, summarize, transpose the dataframe in many ways
# usually the most difficult part is seeing the steps from a to b 

# we want to analyze the otu data by variables in the mapping file
otus<-all_otus.melt
str(otus)

# updated mapping mapBr2.txt
#map<-read.table("mapBr2.txt")#, col.names = c('sample','Elevation','rep','d'))
str(map)

library(tidyr)# tidyr
# reference - http://tidyr.tidyverse.org/

# tidy data
# each variable is in a column
# each observation is a row
# each value is a cell
View(otus)
head(otus) # not tidy - to be tidy, each row = sample and each col = taxa (at desired grouping) and cell = abundance (count)
# this is our goal

# 4-fundamental functions of data tidying - gather, separate, unite, spread
# unite() combines multiple columns into a single column
# separate() splits a single column into multiple columns

#combine taxonomy columns
colnames(otus)
#otus<-unite(otus, col='taxa', kingdom:species, remove=TRUE)
#colnames(otus) #"remove = TRUE" dropped the former columns
View(otus)
# reverse:
#otu_undo<-separate(otus, taxa, into=c('kingdom','phylum','class','order','family','genus','species'), sep=';', fill='right')
#colnames(otu_undo)
#rm(otu_undo)# remove from env

# gather() takes multiple columns, and gathers them into key-value pairs: it makes “wide” data longer
# spread() takes two columns (key & value) and spreads in to multiple columns, it makes “long” data wider

######################################
##DID NOT TRY AND USE ANY SCRIPTS BELOW JULY 2018############

######################################


## creating new variables, working with strings and factors
# string = 'character'
# single or double quotes - no difference

# we need to link the data from the table to to tie to mapping data
# make a new variable that tells the chamberstatus
unique(otus$sample) # prints with "" if character string,
unique(otus$SID) # prints with "" if character string,
unique(otus$taxa) # prints with "" if character string,
otus$sample<-as.factor(otus$sample) # levels if factor. character can be similar to factors, but treated differently in r
unique(otus$sample) # can easily convert between
otus$SID<-as.character(otus$SID) # try as.character() when getting errors - fixes problems that occur with factors
unique(otus$taxa)

# we could use tidyr to separate out the species name:
otu_sp<-separate(otus, taxa, sep='', into=c('v1','name','otu','taxon'))
head(otu_sp) # did it but now we have unwanted columns (v1, otu, taxon)
View(otu_sp)
otu_sp_new<-otu_sp[,-c(2,4,5)] # drop extra columns
head(otu_sp_new)
View(otu_sp_new)
#then use gsub to clean '_report' from name
#otu_sp_new$name<-gsub(pattern='_report',replacement='', otu_sp_new$name)
#unique(otu_sp_new$name) 
#close

#alternatively, we can use 3 'gsub'instead of using separate
?gsub #replaces all pattern matches in a string

# each file is sandwiched by 'dustProject/' , '/otu/otu_taxon.xls', and '_report'
otus$sample<-gsub(pattern='MM', replacement='', otus$sample)#
View(otus)
unique(otus$sample)
#otus$name<-gsub(pattern='/otu/otu_taxon.xls', replacement='', otus$name)
#otus$name<-gsub(pattern='_report', replacement='', otus$name)
#unique(otus$name) #nearly same result BUT different syntax for 1 folder (report 2) isn't as clean

# this is where consistency in naming things can make things easier
# or just due different method -
#levels(otus$file)<-c('M16s','MITS1','TITS1','T16s')
head(otus)

#now want a separate column for the gene region and Month species- currently combined
# mutuate from 'dplyr' can be used for new variables

######################################
# dplyr & reshape2
# set operations, subsetting, filtering, summarizing, joining
# ref - http://dplyr.tidyverse.org/articles/dplyr.html

# remove where count = 0 (not present) to reduce df size
otu.count<-otus%>%filter(count >= 1)
str(otu.count)

# now ready to reshape data, group & summarize by desired groupings
# currently in long format - row for each sample x taxa
# group by taxanomic ID
otu.taxa<-otu.count%>%
  group_by(sample, taxa, SID)%>%
  summarize(sum_otu = sum(count)) #taking a long time
##not sure if denovo ID is speficic within Month/gene groups or universal for all otu ids
# should generate df with the abundance of each taxanomic group for each sample
View(otu.taxa)  

length(unique(otu.taxa$taxa))#different Taxa ids - 372
#rename for each of use or combine with data containing groupings

# join map to otus
# joins - dplyr - http://dplyr.tidyverse.org/articles/two-table.html
# left_join(left, right, by=c(left = right)) OR by = both if varible has same name in both dfs
# right_join(left, right, by=c(left, right))
#otu.guilds<-otu.taxa%>%
# filter(gene=='ITIS')%>%
# left_join(guilds, by=c('ID' = 'OTU.ID'))
#View(otu.guilds)
# not working - i think ened to use the OTU.ID to join
# for now - let's just use the order

##generate otu guild with combined otu table
#grps<-otu.taxa%>%
# filter(gene=='ITIS')%>%
#group_by(taxa)%>%
#summarize(sum=sum(sum_otu), n_otus=length(unique(ID)))

grps<-otu.taxa%>%
  group_by(taxa)%>%
  summarize(sum=sum(sum_otu), n_otus=length(unique(SID)))

length(unique(grps$taxa))# need corresponding 372 otus groups

#otu.order<-otu.count%>%
# separate(taxa, into=c('kingdom','phylum','class','order','family','genus','species'), sep=';')%>%
#dplyr::select(order,sample, count)
#str(otu.order)
#View(otu.order)

otu.class<-otu.count%>%
  separate(taxa, into=c('kingdom','phylum','class','order','family','genus','species'), sep=';')%>%
  dplyr::select(class,sample, count)
str(otu.class)
View(otu.class)

length(unique(otu.class$class)) #63 different classes

#clean up to make sure is accurate
otu.class$class<-tolower(otu.class$class)
otu.class$class<-gsub('c', '', otu.class$class)
otu.class$class<-gsub('_', '', otu.class$class)
length(unique(otu.class$class)) #63

# cast dataframe so there is a column for each order - reshape2
# guide - http://seananderson.ca/2013/10/19/reshape.html
# melt and dcast are opposites
otu.cast<-otu.class%>%
  dcast(sample~class, value.var='count')
#error - aggregation function missing defalting to length
# this occurs when there are multiple rows for each sample grouping
# in this case - each row is an otu/denovo and there may be multiple within each order

# summarize it, then cast
otu.cast<-otu.class%>%
  group_by(sample, class)%>% #use same groupings as waned in the cast
  summarize(sum = sum(count))%>%
  dcast(sample~class)
dim(otu.cast) #14 rows of 64 variables

# now to run multivariate tests, want to combine with mapping data tied to each sample
# in 'otu.df', each 'sample' corresponds to different treatment levels in the 'map' df
# we can use dplyr joins to add those to the corresponding data points

# meta data 'TreatName' matches to 'sample' in otu.cast
unique(map$sample)
unique(otu.cast$sample)
# different data types
map$sample<-as.factor(map$sample)

# use set operations to make sure all levels match
diff<-setdiff(map$sample, otu.cast$sample) # returns observations in x (map.meta) but not in y (otu.df)
diff
#"E15010" "E15011" "E15012"

#levels that are the same in both dfs
same<-intersect(map$sample, otu.cast$sample) #"CF4"    "CF5"    "CF6"    "E100R7" "E100R8" "E100R9"
same
# 6 of the names match, 3 do not - in the otu.df there is an 'R' in the E150 names
unique(otu.cast$sample)
unique(map$sample)

#otu.cast$treat<-ifelse(otu.cast$sample %in% same, otu.cast$sample,#if sample mathes to TreatName in map.meta keep that sample
#gsub(pattern='R',replacement='', otu.cast$sample))# for all other rows, replace the 'R' with nothing('') from sample
#unique(otu.cast$treat) # created a mix of numbers and the 'E150' names without R
# b/c -the data types got mixed up - using gsub creates a character vector but inputs are factors
# converts factors to numbers if not designated as character

#otu.cast$treat<-ifelse(otu.cast$sample %in% same, as.character(otu.cast$sample),#add as.character
#gsub(pattern='R',replacement='', otu.cast$sample))
#unique(otu.cast$treat)

#compare mapping and otu again
is.element(unique(otu.cast$sample), unique(map$sample)) #gives TRUE or FALSE if each level matches
unique(otu.cast$sample)

# join map to otus
# joins - dplyr - http://dplyr.tidyverse.org/articles/two-table.html
# left_join(left, right, by=c(left, right)) OR by = both if varible has same name in both dfs
# right_join(left, right, by=c(left, right))
otu_map<-left_join(otu.cast, map, by=c('sample'='sample'))
# the corresponding columns in these dfs are otus$sample and map$TreatName as show in the 'by'
# the warning that comes up is because the data types are different, but is OK

str(otu_map) #twice as many rows as the otu.cast so somethign was duplicated
View(otu_map)
# created 'Month.x' and 'Month.y' columns
# if there is a shared column name that is not designated in 'by', will avoid duplicate column names by adding '.x' etc
# these are important for the join -
#otu_map2<-left_join(otu.cast, map, by=c('sample'='sample', 'Elevation')) # adds map.meta to end
dim(otu_map)# 36 rows

otu_map3<-left_join(otu.cast, map, by=c('sample'))%>%
  dplyr::select(sample, everything())# reorder columns
View(otu_map3)

write.csv(otu_map,'/Users/maltz/Desktop/RdataDust/otus_by_class.csv', row.names=FALSE)
write.csv(otu_map3,'/Users/maltz/Desktop/RdataDust/otus_by_class3.csv', row.names=FALSE)
write.csv(otu.taxa,'/Users/maltz/Desktop/RdataDust/otus.taxa.csv', row.names=FALSE)
write.csv(otu.class,'/Users/maltz/Desktop/RdataDust/otus.class.csv', row.names=FALSE)

#########################################################
######################################
## calculate dissimilarity/distnace matrics, weighted and unweighted
# community composition differences between dust treatment groups, also Month species and interaction
otu_map<-read.csv('/Users/maltz/Desktop/RdataDust/otus_by_class.csv')

# for community analyses typically need community matrix, maybe dissimilarity matrix, grouping variables
colnames(otu_map)
comm.mat<-otu_map%>%dplyr::select([hlraidbateria]:`zb2`) # select columns for classes
comm.mat[is.na(comm.mat)]<- 0 #replaces 'NA's with 0 for absence

## transforming data prior to dissimilarity
# ecological data commonly has a lot o f 0's, or many rare species and a few super abundant
# transform to account for this
library(vegan)
?decostand # explanation of standardizations under 'details'

#relative composition
comm.rel<-decostand(comm.mat, method='total')
rowSums(comm.rel) #each sample totals to 1 - abundances are now the proportion of total counts
comm.d<-decostand(comm.mat, method='total', MARGIN=2) #defaults to MARGIN = 1
rowSums(comm.d) #not 1
colSums(comm.d) #1 - proportions by column - i.e. how are orderd distributed among samples

# wisconsin double standardization
# common in ecology
# columns are standardized by the max in each column, and then for each sample by the total
# draws out variations in less common taxa
comm.wi<-wisconsin(comm.mat)

# dissimilarity
?vegdist #options for dissimilarity indices under 'method'
comm.dis<-vegdist(log(1+comm.mat), method='bray')
class(comm.dis) #distance matrix

## dissimilarity methods
# bray - includes information about abundances in dissimilarity calculation, non-euclidean
#sum of absolute difference in counts divided by sum of abundances between two samples
# i.e. #different/total#
# differences in methods usually based around how shared 0's are treated, how abundance is considered
comm.jacc<-vegdist(comm.mat, method='jaccard') # treats as 0&1s even if not

## permanova
dust.perm<-adonis(comm.dis~pH+Month*dust, data=otu_map, permutations=10000)
dust.perm

# which variables best explain differences in compositoin?
#rda(community matrix ~ constraining vars)
comm.pca<-rda(comm.mat)
head(summary(comm.pca)) #without any constraining variables this is a pca
plot(comm.pca) # 2 major groups separate out

#add constraining variables - what best explains variation in community composition?
comm.con<-rda(comm.mat~Month+dust+pH+TotalP+OrganicMatter+AvailK, data=otu_map)
#condition - partial out the effects of something
head(summary(comm.con))# produces biplot scores for variables
plot(comm.con)




######################################
## Part 2 - data wrangling 
# By CN
# 02/07/2018

# Preliminaries
#setwd('/Users/maltz/Desktop/RdataOz')

list.files('dustProject')
library(readr)
all_otus_melt <- read_csv("all_otus_melt.csv")
View(all_otus_melt)
library(readxl)
library(dplyr)
library(reshape2)
library(tidyr)
library(vegan)
library(ggplot2)# ggplot resource -http://rpubs.com/collnell/ggplot2

library(tidyverse) # useful packages in 1 - dplyr, ggplot2, tidyr +


######################################
## creating statistical pipelines 
# making workflow and analyses reproducbile, efficient, easy to repurpose

## building pipelines - think of everything in pieces, checking at each step to make sure has desired effect
# then, reconstruct working pieces into short & understandable chunks, clean it up
# often not obvious where steps can be simplified or piped together during first work through, ordering
# if code is not working - put '#' after lines to locate problem 

## remember - '%>%" is called a pipe or a chain
# it chains together a series of commands administered to the same data frame
# common across 'tidyverse' packages and many more
# benefit - use consistent syntax, organized and easy to read, do many things in one chunk, reduces typing & related errors

## cleaning up progress made last week (dust_wrangling.R) 
# streamline data cleaning workflow

## read in data and clean up variable codings
library(data.table) #fread function

otus<-fread('all_otus_melt.csv')%>% # reads in large datasets faster than 'read.csv'
  filter(count >= 1)%>% # drop 0 abundance data
  unite(col='taxa', kingdom:species, remove=TRUE)%>%
  mutate(file = recode(file,"dustProject/Michangensis16s_report/otu/otu_taxon.xls"='M_16s', # manually recoding the file name to reflect Month and gene
                       "dustProject/MichangensisITS1_report 2/otu/otu_taxon.xls" = 'M_ITS1',
                       "dustProject/Tchinensis-ITS1_report/otu/otu_taxon.xls" = 'T_ITS1', 
                       "dustProject/Tchinensis16s_report/otu/otu_taxon.xls" = 'T_16s'))%>%
  separate(col='file', into=c('Month','gene'), sep='_')%>% # separate into separate columns for Month and gene
  group_by(Month, gene, sample, taxa, ID)%>% # set grouping level you want data aggregated at
  summarize(otu_abun = sum(count), otu_rich =length(unique(ID))) #%>% # summarize the abundances at desired level - here is just each unique ID x taxa (so every otu)

## otu abundance vs otu richness
# otu_abun =e total abundance (counting multiple reads of the same otu), otu_rich = # of unique IDs - denovos
# here, otu_rich should be 1 because we are summarizing at the ID level, BUT if done differently, will tell you the richness of otus vs the abundance
# worth considerating - what is more meaningful, total reads or unique otus? are reads skewed by methodology or do they reflect abundance in an accurate way?
# also - how to deal with multiple 'ID's of same 'taxa'? read as a distinct otu, but taxanomic difference unknown - consistent across Month species

head(otus) # in long format, where you can easily bind with alternative grouping data (functional groups)
write.csv(otus, 'dustProject/otus_by_ID.csv', row.names=FALSE)

## make community data matrix for analysis
# summarize data at desired groupings for community data matrix (what do you want your 'species' or columns to be?)
## If you are going to regroup that into another taxonomic split -- species  
# summarize at the order level to look at general patterns in community composition
# for other - can use 'left_join' to bind data, then do these next steps

otu.order<-otus%>%
  filter(gene=='ITS1')%>% # use only one gene
  separate(taxa, into=c('kingdom','phylum','class','order','family','genus','species'), sep=';')%>% #this produces a warning 'too few values' because not all taxa are id'd to species, not an issue
  dplyr::select(order, sample:Month, otu_abun, otu_rich)%>% # drop other taxanomic columns not using
  mutate(order = tolower(as.character(order)))%>%
  mutate(order = gsub('_o', '', order))%>% # clean up order names for consistency
  mutate(order = gsub('_', '', order))%>% 
  group_by(sample, gene, Month, order)%>% # the lowest grouping level here (order) is what you want your columns to be
  summarize(otu_abun= sum(otu_abun), otu_rich = sum(otu_rich))#%>% # summarize at the order level
head(otu.order) # df that tells you the otu abundance and richness for each sample of each Month species at the order level (long format)

# note - different than before I am using the 'mutate' function with 'tolower', 'gsub' and 'recode' within them
# this works b/c mutate makes new variables, and those other commands work on vectors, not dataframes
# makes it easier to pipe together for flow and don't need to write 'otu.order$' in front of everything

unique(otu.order$order)# use to verify order names - duplicates? misspellings? NA are those that didn't get id'd to order

## cast into community dataframe (samples as rows, species as columns), & add grouping variables (WIDE DATA)
# here need to decide whether total otu abundance or the richness of reads matters more - use to cast
otu.cast<-otu.order%>%
  dcast(sample+Month+gene~order, value.var='otu_abun')#%>% # for each sample (lowest level within Month, gene, sample), otu abundance of each order (as columns)
# replace NA with 0 in abundances
otu.cast[is.na(otu.cast)]<- 0 
head(otu.cast)

## read in mapping data
# remember - not all sample and treatnames matched due to an R in some names
map.meta<-read_excel("dustProject/dust_map4_wMetadata.xlsx", sheet=2)%>%
  dplyr::select(-contains('Sequence'), -Variable, -Description, -`#SampleID`)%>% # drop cols not needed
  mutate(treat = recode(TreatName, "E15010"='E150R10', "E15011"="E150R11", "E15012"="E150R12", .default = TreatName)) # default gives variable value if not specifically named in 'recode'
View(map.meta)
unique(map.meta$treat)
colnames(map.meta)
write.csv(map.meta, 'dustProject/dust_metadata.csv', row.names=FALSE)# rewite so don't need to recode variables again

# combine
otu_map<-left_join(otu.cast, map.meta, by=c('sample'='treat', 'Month'))%>%
  dplyr::select(sample, gene, Month, dust:SoilDryWeight, everything())
View(otu_map)

##This only has Month T
write.csv(otu_map, 'dustProject/otus_by_order.csv', row.names=FALSE)
View










################################################################################################
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ################################################################################################  
## goals - read in mapping file, OTU tables, and distance matrices (weighted or unweighted)  
# assess community composition differences between chamber treatment groups
# work with OTU table which used numeric values, (as opposed to presence absence) or percentages; 
# expand stats using distance matrices. 
# PERMANOVA,exploring other community analyses, RDA
## preliminaries 
library(dplyr)
## data
# Two chamber statuses, Ch and N
# using files in RdataDust folder (contianing dust)
setwd('/Users/maltz/Desktop/RdataDust')#downlaoded location of files
list.files()#show filenames in current working directory
list.dirs()#everything in folder 

# read in mapping file
list.files()
# column A in the mapping file, ‘sample’ hmnas of the samples as given on the OTU table,
map<-read_excel("SierraMap6.xlsx", col.names =c('SampleID','Year','Month','Site','SiteRep','Elevation','DateCode','Description'))
map1<-read_excel("SierraMap6.xlsx", col_names = c('SampleID','BarcodeSequence','LinkerPrimer','Year','Month','SiteCode','RepNum','SiteRep','Site','Elevation','DateCode','DescName','Description'))
str(map1)

map2<-read_excel("SierraMap6.xlsx", col_names = c('SampleID','BarcodeSequence','LinkerPrimer','Year','Month','SiteCode','RepNum','SiteRep','Site','Elevation','DateCode','DescName','Description'))
str(map2)

map3<-read.table("SierraMap6.txt", header =TRUE, col.names =c('SampleID','Year','Month','Site','SiteRep','Elevation','DateCode','Description','d'))
#map<-read.table("mapSP9.txt")
View(map3)

# read in otu table
list.files()
## reading in file & common error troubleshooting 
library(readxl) #read xlsx files into r on mac computer
# guide: https://github.com/tidyverse/readxl
# if you only have a few files, it is easy enough to read them in individually
# read in single file of '.xls' format
dust_otu<-read_excel("dust_otu.xlsx")
#count.fields(dust_otu)
#Botu<-read_excel("dust_otu.xlsx")
# not working for some reason - try as table
#b_otu<-read.table("dust_otu.txt")
# Error:line 43 did not have 22 elements
# usually means that column headers don't match up with the number of columns

# how many columns are in each row? 
count.fields("dust_otu.txt") # number of cols in each row are not equal
# this is b/c delimination is not clear - in this case b/c cells are tabbed &
# R is reading in the 'taxonomy' column to have different lengths based on resolution

?read.table
# fix - use 'fill' to fill in missing cells so all rows are equal
#b_otu<-read.table("dust_otu.txt", fill=TRUE)
#View(b_otu) # but the first line is the column header
dust_otu.txt<-read.table("dust_otu.txt", fill=TRUE, header=TRUE)
# when rows are equal, this designates the column names from the first row, but not here
View(dust_otu.txt)
# instead, skip line 1 and assign own names
dust_otu<-read.table("dust_otu.txt", fill=TRUE, skip=1, 
                        col.names=c('SID','MMCH5','MMCH5R','MMCH1R','MMCH1','MMN4','MMN3R','MMN4R','MMCH3R','MMCH4','MMN3','MMN2','MMCH3','MMN2R','MMCH4R',
                                    'kingdom','phylum','class','order','family','genus','species'))
View(dust_otu)# great...but what if it is too tedious to type out all the column names?

View(taxa_list)

all_otus<-dust_otu
View(all_otus)# now they are all in a single file

library(dplyr)
write.csv(all_otus, '/Users/maltz/Desktop/RdataDust/all_otus5.csv', row.names=FALSE)
list.files()
# now this last part (binding the dataframes) will only work if the column headers are exactly the same in each file
# so if there were different numbers of samples between the files...

# melt dataframe so there is a column for each sample name
# basically reduce a subset of columns into 2 columns (1 with the former column name, 1 with the values)
library(reshape2)
#View(file2)
#View(all_otus)
#all_otus.melt<-melt(all_otus, id.vars=c('SID',taxa_list))
#View(all_otus.melt) 
#####################################################
otu.class<-otus%>%
  separate(taxa, into=c('kingdom','phylum','class','order','family','genus','species'), sep=';')%>% #this produces a warning 'too few values' because not all taxa are id'd to species, not an issue
  dplyr::select(class, sample, otu_abun, otu_rich)%>% # drop other taxanomic columns not using
  mutate(class = tolower(as.character(class)))%>%
  mutate(class = gsub('_c', '', class))%>% # clean up class names for consistency
  mutate(class = gsub('_', '', class))%>% 
  group_by(sample, class)%>% # the lowest grouping level here (class) is what you want your columns to be
  summarize(otu_abun= sum(otu_abun), otu_rich = sum(otu_rich))#%>% # summarize at the class level
head(otu.class) # df that tells you the otu abundance and richness for each sample of each Month species at the class level (long format)
write.csv(otu.class, '/Users/maltz/Desktop/RdataDust/otu_class.csv', row.names=FALSE)
# note - different than before I am using the 'mutate' function with 'tolower', 'gsub' and 'recode' within them
# this works b/c mutate makes new variables, and those other commands work on vectors, not dataframes
# makes it easier to pipe together for flow and don't need to write 'otu.class$' in front of everything

unique(otu.class$class)# use to verify class names - duplicates? misspellings? NA are those that didn't get id'd to class
########################################
# include syntax to name new columns
#all_otus.melt<-melt(all_otus, id.vars=c('SID',taxa_list), value.name='count', variable.name='sample')
#str(all_otus.melt)
# to consider another time - 
# how to deal with otus that are same taxomony but different otu?
# are the otu id's 'denovo0' the same otu across files? or these are unique to each file?

## add 'melt' to function
#header_melt_fun<-function(x){
# taxa_list<-c('kingdom','phylum','class','order','family','genus','species') #taxonomy column headers
#temp.df<-read.table(otu_path[1], fill=TRUE, stringsAsFactors = FALSE)
#head_list<-as.character(temp.df[1, 1:(length(temp.df)-7)]) #dataframe headers
#out.df<-read.table(x, fill=TRUE, skip=1,col.names=c(head_list,taxa_list)) # read in and assign headers
#out.df$file<-paste(x) # create additional column that designates source file
#out.df.melt<-melt(out.df, id.vars=c('ID','file', taxa_list), value.name='count', variable.name='sample')
#return(out.df.melt)
#}

#all_otus_melt2<-bind_rows(lapply(otu_path, header_melt_fun))
#View(all_otus.melt) #everything's here

#write.csv(all_otus_melt, '/Users/maltz/Desktop/RdataDust/all_otus_melt4.csv', row.names=FALSE)

######################################
## wranging skills
# manipulating data with dplyr, tidyr, and reshape2, working with characeter strings
# familiarity with different functions, splitting and combinging existing variables

setwd('/Users/maltz/Desktop/RdataDust')# downlaoded location of files
list.files()

library(dplyr)
library(reshape2)

# from otus created above, we can regroup, summarize, transpose the dataframe in many ways
# usually the most difficult part is seeing the steps from a to b 

# we want to analyze the otu data by variables in the mapping file
otus<-all_otus.melt
str(otus)

# updated mapping mapBr2.txt
#map<-read.table("mapBr2.txt")#, col.names = c('sample','Elevation','rep','d'))
str(map)

library(tidyr)# tidyr
# reference - http://tidyr.tidyverse.org/

# tidy data
# each variable is in a column
# each observation is a row
# each value is a cell
#View(otus)
#head(otus) # not tidy - to be tidy, each row = sample and each col = taxa (at desired grouping) and cell = abundance (count)
# this is our goal

# 4-fundamental functions of data tidying - gather, separate, unite, spread
# unite() combines multiple columns into a single column
# separate() splits a single column into multiple columns

#combine taxonomy columns
#colnames(otus)
#otus<-unite(otus, col='taxa', kingdom:species, remove=TRUE)
#colnames(otus) #"remove = TRUE" dropped the former columns
#View(otus)
# reverse:
#otu_undo<-separate(otus, taxa, into=c('kingdom','phylum','class','order','family','genus','species'), sep=';', fill='right')
#colnames(otu_undo)
#rm(otu_undo)# remove from env

# gather() takes multiple columns, and gathers them into key-value pairs: it makes “wide” data longer
# spread() takes two columns (key & value) and spreads in to multiple columns, it makes “long” data wider

######################################
## creating new variables, working with strings and factors
# string = 'character'
# single or double quotes - no difference

# we need to link the data from the table to to tie to mapping data
# make a new variable that tells the chamberstatus
#unique(otus$sample) # prints with "" if character string,
#unique(otus$SID) # prints with "" if character string,
#unique(otus$taxa) # prints with "" if character string,
#otus$sample<-as.factor(otus$sample) # levels if factor. character can be similar to factors, but treated differently in r
#unique(otus$sample) # can easily convert between
#otus$SID<-as.character(otus$SID) # try as.character() when getting errors - fixes problems that occur with factors
#unique(otus$taxa)

# we could use tidyr to separate out the species name:
#otu_sp<-separate(otus, taxa, sep='', into=c('v1','name','otu','taxon'))
#head(otu_sp) # did it but now we have unwanted columns (v1, otu, taxon)
#View(otu_sp)
#otu_sp_new<-otu_sp[,-c(2,4,5)] # drop extra columns
#head(otu_sp_new)
#View(otu_sp_new)
#then use gsub to clean '_report' from name
#otu_sp_new$name<-gsub(pattern='_report',replacement='', otu_sp_new$name)
#unique(otu_sp_new$name) # "Michangensis16s"    "MichangensisITS1 2" "Tchinensis-ITS1"    "Tchinensis16s"
#close

#alternatively, we can use 3 'gsub'instead of using separate
?gsub #replaces all pattern matches in a string

# each file is sandwiched by 'dustProject/' , '/otu/otu_taxon.xls', and '_report'
otus$sample<-gsub(pattern='MM', replacement='', otus$sample)#
View(otus)
unique(otus$sample)
#otus$name<-gsub(pattern='/otu/otu_taxon.xls', replacement='', otus$name)
#otus$name<-gsub(pattern='_report', replacement='', otus$name)
#unique(otus$name) #nearly same result BUT different syntax for 1 folder (report 2) isn't as clean

# this is where consistency in naming things can make things easier
# or just due different method -
#levels(otus$file)<-c('M16s','MITS1','TITS1','T16s')
head(otus)

#now want a separate column for the gene region and Month species- currently combined
# mutuate from 'dplyr' can be used for new variables

######################################
# dplyr & reshape2
# set operations, subsetting, filtering, summarizing, joining
# ref - http://dplyr.tidyverse.org/articles/dplyr.html

# remove where count = 0 (not present) to reduce df size
otu.count<-otus%>%filter(count >= 1)
str(otu.count)

# now ready to reshape data, group & summarize by desired groupings
# currently in long format - row for each sample x taxa
# group by taxanomic ID
otu.taxa<-otu.count%>%
  group_by(sample, taxa, SID)%>%
  summarize(sum_otu = sum(count)) #taking a long time
##not sure if denovo ID is speficic within Month/gene groups or universal for all otu ids
# should generate df with the abundance of each taxanomic group for each sample
View(otu.taxa)  

length(unique(otu.taxa$taxa))#different Taxa ids - 372
#rename for each of use or combine with data containing groupings

# join map to otus
# joins - dplyr - http://dplyr.tidyverse.org/articles/two-table.html
# left_join(left, right, by=c(left = right)) OR by = both if varible has same name in both dfs
# right_join(left, right, by=c(left, right))
#otu.guilds<-otu.taxa%>%
# filter(gene=='ITIS')%>%
# left_join(guilds, by=c('ID' = 'OTU.ID'))
#View(otu.guilds)
# not working - i think ened to use the OTU.ID to join
# for now - let's just use the order

##generate otu guild with combined otu table
#grps<-otu.taxa%>%
# filter(gene=='ITIS')%>%
#group_by(taxa)%>%
#summarize(sum=sum(sum_otu), n_otus=length(unique(ID)))

grps<-otu.taxa%>%
  group_by(taxa)%>%
  summarize(sum=sum(sum_otu), n_otus=length(unique(SID)))

length(unique(grps$taxa))# need corresponding 372 otus groups

#otu.order<-otu.count%>%
# separate(taxa, into=c('kingdom','phylum','class','order','family','genus','species'), sep=';')%>%
#dplyr::select(order,sample, count)
#str(otu.order)
#View(otu.order)

otu.class<-otu.count%>%
  separate(taxa, into=c('kingdom','phylum','class','order','family','genus','species'), sep=';')%>%
  dplyr::select(class,sample, count)
str(otu.class)
View(otu.class)

length(unique(otu.class$class)) #63 different classes

#clean up to make sure is accurate
otu.class$class<-tolower(otu.class$class)
otu.class$class<-gsub('c', '', otu.class$class)
otu.class$class<-gsub('_', '', otu.class$class)
length(unique(otu.class$class)) #63

# cast dataframe so there is a column for each order - reshape2
# guide - http://seananderson.ca/2013/10/19/reshape.html
# melt and dcast are opposites
otu.cast<-otu.class%>%
  dcast(sample~class, value.var='count')
#error - aggregation function missing defalting to length
# this occurs when there are multiple rows for each sample grouping
# in this case - each row is an otu/denovo and there may be multiple within each order

# summarize it, then cast
otu.cast<-otu.class%>%
  group_by(sample, class)%>% #use same groupings as waned in the cast
  summarize(sum = sum(count))%>%
  dcast(sample~class)
dim(otu.cast) #14 rows of 64 variables

# now to run multivariate tests, want to combine with mapping data tied to each sample
# in 'otu.df', each 'sample' corresponds to different treatment levels in the 'map' df
# we can use dplyr joins to add those to the corresponding data points

# meta data 'TreatName' matches to 'sample' in otu.cast
unique(map$sample)
unique(otu.cast$sample)
# different data types
map$sample<-as.factor(map$sample)

# use set operations to make sure all levels match
diff<-setdiff(map$sample, otu.cast$sample) # returns observations in x (map.meta) but not in y (otu.df)
diff
#"E15010" "E15011" "E15012"

#levels that are the same in both dfs
same<-intersect(map$sample, otu.cast$sample) #"CF4"    "CF5"    "CF6"    "E100R7" "E100R8" "E100R9"
same
# 6 of the names match, 3 do not - in the otu.df there is an 'R' in the E150 names
unique(otu.cast$sample)
unique(map$sample)

#otu.cast$treat<-ifelse(otu.cast$sample %in% same, otu.cast$sample,#if sample mathes to TreatName in map.meta keep that sample
#gsub(pattern='R',replacement='', otu.cast$sample))# for all other rows, replace the 'R' with nothing('') from sample
#unique(otu.cast$treat) # created a mix of numbers and the 'E150' names without R
# b/c -the data types got mixed up - using gsub creates a character vector but inputs are factors
# converts factors to numbers if not designated as character

#otu.cast$treat<-ifelse(otu.cast$sample %in% same, as.character(otu.cast$sample),#add as.character
#gsub(pattern='R',replacement='', otu.cast$sample))
#unique(otu.cast$treat)

#compare mapping and otu again
is.element(unique(otu.cast$sample), unique(map$sample)) #gives TRUE or FALSE if each level matches
unique(otu.cast$sample)

# join map to otus
# joins - dplyr - http://dplyr.tidyverse.org/articles/two-table.html
# left_join(left, right, by=c(left, right)) OR by = both if varible has same name in both dfs
# right_join(left, right, by=c(left, right))
otu_map<-left_join(otu.cast, map, by=c('sample'='sample'))
# the corresponding columns in these dfs are otus$sample and map$TreatName as show in the 'by'
# the warning that comes up is because the data types are different, but is OK

str(otu_map) #twice as many rows as the otu.cast so somethign was duplicated
View(otu_map)
# created 'Month.x' and 'Month.y' columns
# if there is a shared column name that is not designated in 'by', will avoid duplicate column names by adding '.x' etc
# these are important for the join -
#otu_map2<-left_join(otu.cast, map, by=c('sample'='sample', 'Elevation')) # adds map.meta to end
dim(otu_map)# 36 rows

otu_map3<-left_join(otu.cast, map, by=c('sample'))%>%
  dplyr::select(sample, everything())# reorder columns
View(otu_map3)

write.csv(otu_map,'/Users/maltz/Desktop/RdataDust/otus_by_class.csv', row.names=FALSE)
write.csv(otu_map3,'/Users/maltz/Desktop/RdataDust/otus_by_class3.csv', row.names=FALSE)
write.csv(otu.taxa,'/Users/maltz/Desktop/RdataDust/otus.taxa.csv', row.names=FALSE)
write.csv(otu.class,'/Users/maltz/Desktop/RdataDust/otus.class.csv', row.names=FALSE)

#########################################################
######################################
## calculate dissimilarity/distnace matrics, weighted and unweighted
# community composition differences between dust treatment groups, also Month species and interaction
otu_map<-read.csv('/Users/maltz/Desktop/RdataDust/otus_by_class.csv')

# for community analyses typically need community matrix, maybe dissimilarity matrix, grouping variables
colnames(otu_map)
comm.mat<-otu_map%>%dplyr::select([hlraidbateria]:`zb2`) # select columns for classes
comm.mat[is.na(comm.mat)]<- 0 #replaces 'NA's with 0 for absence

## transforming data prior to dissimilarity
# ecological data commonly has a lot o f 0's, or many rare species and a few super abundant
# transform to account for this
library(vegan)
?decostand # explanation of standardizations under 'details'

#relative composition
comm.rel<-decostand(comm.mat, method='total')
rowSums(comm.rel) #each sample totals to 1 - abundances are now the proportion of total counts
comm.d<-decostand(comm.mat, method='total', MARGIN=2) #defaults to MARGIN = 1
rowSums(comm.d) #not 1
colSums(comm.d) #1 - proportions by column - i.e. how are orderd distributed among samples

# wisconsin double standardization
# common in ecology
# columns are standardized by the max in each column, and then for each sample by the total
# draws out variations in less common taxa
comm.wi<-wisconsin(comm.mat)

# dissimilarity
?vegdist #options for dissimilarity indices under 'method'
comm.dis<-vegdist(log(1+comm.mat), method='bray')
class(comm.dis) #distance matrix

## dissimilarity methods
# bray - includes information about abundances in dissimilarity calculation, non-euclidean
#sum of absolute difference in counts divided by sum of abundances between two samples
# i.e. #different/total#
# differences in methods usually based around how shared 0's are treated, how abundance is considered
comm.jacc<-vegdist(comm.mat, method='jaccard') # treats as 0&1s even if not

## permanova
dust.perm<-adonis(comm.dis~pH+Month*dust, data=otu_map, permutations=10000)
dust.perm

# which variables best explain differences in compositoin?
#rda(community matrix ~ constraining vars)
comm.pca<-rda(comm.mat)
head(summary(comm.pca)) #without any constraining variables this is a pca
plot(comm.pca) # 2 major groups separate out

#add constraining variables - what best explains variation in community composition?
comm.con<-rda(comm.mat~Month+dust+pH+TotalP+OrganicMatter+AvailK, data=otu_map)
#condition - partial out the effects of something
head(summary(comm.con))# produces biplot scores for variables
plot(comm.con)




######################################
## Part 2 - data wrangling 
# By CN
# 02/07/2018

# Preliminaries
#setwd('/Users/maltz/Desktop/RdataOz')

list.files('dustProject')
library(readr)
all_otus_melt <- read_csv("all_otus_melt.csv")
View(all_otus_melt)
library(readxl)
library(dplyr)
library(reshape2)
library(tidyr)
library(vegan)
library(ggplot2)# ggplot resource -http://rpubs.com/collnell/ggplot2

library(tidyverse) # useful packages in 1 - dplyr, ggplot2, tidyr +


######################################
## creating statistical pipelines 
# making workflow and analyses reproducbile, efficient, easy to repurpose

## building pipelines - think of everything in pieces, checking at each step to make sure has desired effect
# then, reconstruct working pieces into short & understandable chunks, clean it up
# often not obvious where steps can be simplified or piped together during first work through, ordering
# if code is not working - put '#' after lines to locate problem 

## remember - '%>%" is called a pipe or a chain
# it chains together a series of commands administered to the same data frame
# common across 'tidyverse' packages and many more
# benefit - use consistent syntax, organized and easy to read, do many things in one chunk, reduces typing & related errors

## cleaning up progress made last week (dust_wrangling.R) 
# streamline data cleaning workflow

## read in data and clean up variable codings
library(data.table) #fread function

otus<-fread('all_otus_melt.csv')%>% # reads in large datasets faster than 'read.csv'
  filter(count >= 1)%>% # drop 0 abundance data
  unite(col='taxa', kingdom:species, remove=TRUE)%>%
  mutate(file = recode(file,"dustProject/Michangensis16s_report/otu/otu_taxon.xls"='M_16s', # manually recoding the file name to reflect Month and gene
                       "dustProject/MichangensisITS1_report 2/otu/otu_taxon.xls" = 'M_ITS1',
                       "dustProject/Tchinensis-ITS1_report/otu/otu_taxon.xls" = 'T_ITS1', 
                       "dustProject/Tchinensis16s_report/otu/otu_taxon.xls" = 'T_16s'))%>%
  separate(col='file', into=c('Month','gene'), sep='_')%>% # separate into separate columns for Month and gene
  group_by(Month, gene, sample, taxa, ID)%>% # set grouping level you want data aggregated at
  summarize(otu_abun = sum(count), otu_rich =length(unique(ID))) #%>% # summarize the abundances at desired level - here is just each unique ID x taxa (so every otu)

## otu abundance vs otu richness
# otu_abun =e total abundance (counting multiple reads of the same otu), otu_rich = # of unique IDs - denovos
# here, otu_rich should be 1 because we are summarizing at the ID level, BUT if done differently, will tell you the richness of otus vs the abundance
# worth considerating - what is more meaningful, total reads or unique otus? are reads skewed by methodology or do they reflect abundance in an accurate way?
# also - how to deal with multiple 'ID's of same 'taxa'? read as a distinct otu, but taxanomic difference unknown - consistent across Month species

head(otus) # in long format, where you can easily bind with alternative grouping data (functional groups)
write.csv(otus, 'dustProject/otus_by_ID.csv', row.names=FALSE)

## make community data matrix for analysis
# summarize data at desired groupings for community data matrix (what do you want your 'species' or columns to be?)
## If you are going to regroup that into another taxonomic split -- species  
# summarize at the order level to look at general patterns in community composition
# for other - can use 'left_join' to bind data, then do these next steps

otu.order<-otus%>%
  filter(gene=='ITS1')%>% # use only one gene
  separate(taxa, into=c('kingdom','phylum','class','order','family','genus','species'), sep=';')%>% #this produces a warning 'too few values' because not all taxa are id'd to species, not an issue
  dplyr::select(order, sample:Month, otu_abun, otu_rich)%>% # drop other taxanomic columns not using
  mutate(order = tolower(as.character(order)))%>%
  mutate(order = gsub('_o', '', order))%>% # clean up order names for consistency
  mutate(order = gsub('_', '', order))%>% 
  group_by(sample, gene, Month, order)%>% # the lowest grouping level here (order) is what you want your columns to be
  summarize(otu_abun= sum(otu_abun), otu_rich = sum(otu_rich))#%>% # summarize at the order level
head(otu.order) # df that tells you the otu abundance and richness for each sample of each Month species at the order level (long format)

# note - different than before I am using the 'mutate' function with 'tolower', 'gsub' and 'recode' within them
# this works b/c mutate makes new variables, and those other commands work on vectors, not dataframes
# makes it easier to pipe together for flow and don't need to write 'otu.order$' in front of everything

unique(otu.order$order)# use to verify order names - duplicates? misspellings? NA are those that didn't get id'd to order

## cast into community dataframe (samples as rows, species as columns), & add grouping variables (WIDE DATA)
# here need to decide whether total otu abundance or the richness of reads matters more - use to cast
otu.cast<-otu.order%>%
  dcast(sample+Month+gene~order, value.var='otu_abun')#%>% # for each sample (lowest level within Month, gene, sample), otu abundance of each order (as columns)
# replace NA with 0 in abundances
otu.cast[is.na(otu.cast)]<- 0 
head(otu.cast)

## read in mapping data
# remember - not all sample and treatnames matched due to an R in some names
map.meta<-read_excel("dustProject/dust_map4_wMetadata.xlsx", sheet=2)%>%
  dplyr::select(-contains('Sequence'), -Variable, -Description, -`#SampleID`)%>% # drop cols not needed
  mutate(treat = recode(TreatName, "E15010"='E150R10', "E15011"="E150R11", "E15012"="E150R12", .default = TreatName)) # default gives variable value if not specifically named in 'recode'
View(map.meta)
unique(map.meta$treat)
colnames(map.meta)
write.csv(map.meta, 'dustProject/dust_metadata.csv', row.names=FALSE)# rewite so don't need to recode variables again

# combine
otu_map<-left_join(otu.cast, map.meta, by=c('sample'='treat', 'Month'))%>%
  dplyr::select(sample, gene, Month, dust:SoilDryWeight, everything())
View(otu_map)

##This only has Month T
write.csv(otu_map, '/Users/maltz/Desktop/RdataDust'/otus_by_class.csv', row.names=FALSE)
View
