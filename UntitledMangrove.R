####################################
## dust data
## July 04 2018 
## CN for MM
####################################
## goals - read in mapping file, OTU tables, and distance matrices (weighted or unweighted)  
# assess community composition differences between Dust treatment groups (possibly also by Tree species and the interaction), 
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
library(readxl)
View(map)
#map<-map
map1<-read_excel("SierraMap6.xlsx", skip = 1, col_names = c('SampleID','BarcodeSequence','LinkerPrimer','Year','Month','SiteCode','RepNum','SiteRep','Site','Elevation','DateCode','DescName','Description'))

#'#SampleID	BarcodeSequence	LinkerPrimerSequence	Year	Month	SiteCode	RepNum	SiteRep	Site	Elevation	DateCode	DescName	Description
#'
View(map1)

# read in otu tables  for each gene group, and the functional group from FUNGuild
# My Drive > RdataDust > filtered_table_w_metadata.txt
list.files()
#head.list<-read.table("headlist.txt", fill=TRUE, stringsAsFactors = FALSE)
# list of the folders that contain report files
filt<-list.files(pattern ='filtered', full.names=TRUE)
#filt
#list.files(filt) #list the files in those folders
# can see there are same files within each folder
#list.files(filt, pattern='table')

# generate list of the file path to each otu table
otu_path<-paste0(, filtered_table_w_metadata.txt)
otu_path # should be the 4 otu file paths

## reading in files & common error troubleshooting 
# read in otu files
library(readxl) #read xlsx files into r on mac computer
# guide: https://github.com/tidyverse/readxl

# if you only have a few files, it is easy enough to read them in individually

# read in single file of '.xls' format
d16s<-read_excel("filtered_table_w_metadata.xlsx")
View(d16s)
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
# when they are same type of data - i.e. otu tables for each tree & gene region are separate files

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
ldf<-lapply(otu_path, header_fun)
str(ldf)# this produces a list of dataframes
length(ldf) # length is 4 because there were 4 files imported, and a dataframe made for each file
#View(ldf[[2]]) # view second df in list

library(dplyr)
#all_otus<-bind_rows(ldf)##combine list of dataframes
View(dust16s)
all_otus<-dust16s    
View(all_otus)# now they are all in a single file - same result as before

#save
write.csv(all_otus, 'all_otus_16SDust_20180718.csv', row.names=FALSE)
write.csv(all_otus, 'all_otusBD.csv', row.names=FALSE)

# now this last part (binding the dataframes) will only work if the column headers are exactly the same in each file
# so if there were different numbers of samples between the files...

# melt dataframe so there is a column for each sample name
# basically reduce a subset of columns into 2 columns (1 with the former column name, 1 with the values)
library(reshape2)
View(all_otus)

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
all_otus_melt<-bind_rows(lapply(all_otus, header_melt_fun))
View(all_otus_melt) #everything's here

write.csv(all_otus_melt, 'all_otus_melt20180718.csv', row.names=FALSE)



####################################
## Dust data
## Jan 22 2018 
## CN for MM

####################################
## goals - read in mapping file, OTU tables, and distance matrices (weighted or unweighted)  
# assess community composition differences between Dust treatment groups (possibly also by Tree species and the interaction), 
# work with denovo OTU table which used numeric values, (as opposed to presence absence) or percentages; 

# expand stats using distance matrices. 
# PERMANOVA,exploring other community analyses, RDA

## preliminaries 
library(dplyr)

## data
# Two tree species: T (T. chinensis) and M (M. ichangensis)
# Three levels of Dust: 0, 100, 150
# Two gene regions: 16S and ITS

# using files in AcidRainEctoTree_manuscript.materials folder (contianing RdataDust)
setwd('/Users/maltz/Desktop/RdataOz')#downlaoded location of files
# it is possible to get these directly from google drive into R without downloading...can show you this if that is useful for you

####################################
## browsing computer files from in R, reading in, trouleshooting
# familiarize with file paths & types, searching, pitfalls of data import, indexing, lists, 

list.files()#show filenames in current working directory
list.files('RdataDust') #files in 'Dust' folder in the working directory  
list.files('RdataDust')[3] # first file in that folder (for me it's"")
list.files('RdataDust')[1] # first file in that folder (for me it's"16S")
list.files('RdataDust/DustITS1_report 2') # and so on
list.files('RdataDust/DustITS1_report 2', pattern = '.txt') #all txt files

list.dirs('RdataDust')#everything in folder 

# reading in similar file types from different folders
# first - get unique pattern for files
# make a list of all the files
# read them all in

# read in mapping file
list.files('RdataDust')
# From My Drive > Dust Project > Dust_map3.txt
# column H in the mapping file, ‘TreatName’ hmnas of the samples as given on the OTU table in the otu folder,
#map<-read.table("RdataDust/Dust_map3.txt", col.names = c('SampleID','Dust','level','control','tree','TreatName','Description','d'))
map<-read_excel("RdataDust/Dust_map6.xls")#, col.names = c('SampleID','Dust','level','control','tree','TreatName','Variable','AOT40','OrganicMatter','pH','NH4N','NO3N','AmNit','AvailK','AvailP','TotalP','SoilMoisture','GWP','Description','d'))
str(map)
View(map)

# read in otu tables  for each tree species x gene group
# from respective drive folders (Dust16s_report, DustITS_report, , and DustITS_report)
# My Drive > RdataDust > Dust16s_report > otu > otu_taxon.xls
list.files('RdataDust')

# list of the folders that contain report files
folders<-list.files('RdataDust', pattern ='report', full.names=TRUE)
folders
list.files(folders) #list the files in those folders
# can see there are same files within each folder
list.files(folders, pattern='otu')

# generate list of the file path to each otu table
otu_path<-paste0(folders, '/otu/otu_taxon.xls')
otu_path # should be the 4 otu file paths

## reading in files & common error troubleshooting 
# read in 4 otu files
library(readxl) #read xlsx files into r on mac computer
# guide: https://github.com/tidyverse/readxl

# if you only have a few files, it is easy enough to read them in individually

# read in single file of '.xls' format
du16s<-read_excel("RdataDust/Michiganensis16s_report/otu/otu_taxon.xls")
# not working for some reason - try as table
du16s<-read.table("RdataDust/Michiganensis16s_report/otu/otu_taxon.xls")
# Error:line 1 did not have 17 elements
# usually means that column headers don't match up with the number of columns

# how many columns are in each row? 
count.fields("RdataDust/Michiganensis16s_report/otu/otu_taxon.xls") # number of cols in each row are not equal
# this is b/c delimination is not clear - in this case b/c cells are tabbed &
# R is reading in the 'taxonomy' column to have different lengths based on resolution

?read.table
# fix - use 'fill' to fill in missing cells so all rows are equal ##
#fill	
#logical. If TRUE then in case the rows have unequal length, blank fields are implicitly added. See ‘Details’.
du16s<-read.table("RdataDust/Michiganensis16s_report/otu/otu_taxon.xls", fill=TRUE)
View(du16s) # but the first line is the column header
du16s<-read.table("RdataDust/Michiganensis16s_report/otu/otu_taxon.xls", fill=TRUE, header=TRUE)
# when rows are equal, this designates the column names from the first row, but not here

# instead, skip line 1 and assign own names
du16s<-read.table("RdataDust/Michiganensis16s_report/otu/otu_taxon.xls", fill=TRUE, skip=1, 
                  col.names=c('ID','CF4','CF5','CF6','E100R7','E100R8','E100R9','E150R10','E150R11','E150R12',
                              'kingdom','phylum','class','order','family','genus','species'))
View(du16s)# great...but what if it is too tedious to type out all the column names?
##Go back and read in all the OTU tables from each of the three folders
# can save the column names from the file itself
heads<-read.table("RdataDust/Michiganensis16s_report/otu/otu_taxon.xls", fill=TRUE, stringsAsFactors = FALSE)
# using brackets:
# df = [rows,columns]
# can subset any data frame by indicating the rows and cols wanted 
heads[1,] # first row containing the column names
heads[1:9,] # rows 1-9
heads[c(1,2,4,8,9:13),] # rows 1, 2, 4, 8, 9, 10, 11, 12, 13
heads[1:5,1:5]
heads[1:5, -9] #all columns except the 9th
heads[1,1:11] # the given headers are in cols 1:11
heads[1,-c(12:17)] #or all cols but 12:17  - same
head.list<-as.character(heads[1,1:10]) #save row and convert to character string
head.list# using only 1:10 above because we will replace 'taxonomy' with 'kingdom'

# assign the column names to be the given headers + taxonomy 
taxa_list<-c('kingdom','phylum','class','order','family','genus','species') #taxonomy column headers
du16s<-read.table("RdataDust/Michiganensis16s_report/otu/otu_taxon.xls", fill=TRUE, skip=1, 
                  col.names=c(head.list, taxa_list))
View(du16s)
# especially helpful to automate work pipelines and combine multiple data sheets w/ different column names

####################################
## functions, reshaping and combining dataframes, 'apply', automating file imports
# writing a function to read in multiple files, rename columns, and combine 
# when they are same type of data - i.e. otu tables for each tree & gene region are separate files

## goal - read in all 4 otu tables and combine into single df for multivariate tests
otu_path
otu_path[2] # select the 2nd element - when it's a list, here is only 1 dimension in brackets 

# a function to read in otu files and add taxonomy as last 7 columns
header_fun<-function(x){
  taxa_list<-c('kingdom','phylum','class','order','family','genus','species') #taxonomy column headers
  temp.df<-read.table(x, fill=TRUE, stringsAsFactors = FALSE)
  head_list<-as.character(temp.df[1, 1:(length(temp.df)-7)]) #dataframe headers
  out.df<-read.table(x, fill=TRUE, skip=1,col.names=c(head_list,taxa_list)) # read in and assign headers
  out.df$file<-paste(x) # create additional column that designates source file
  return(out.df)
}

#now try this function for a single file path (x)
file1<-header_fun([all_otus])
View(file1)
# other files
file2<-header_fun(otu_path[2])
file3<-header_fun(otu_path[3])
file4<-header_fun(otu_path[4])
#stacks the dataframes on top of each other and combines that and makes it a continuous data frame (##cbind, rows that are the same it will combine them)
all_otus<-rbind(file1, file2, file3, file4)
View(all_otus)# now they are all in a single file

# read in all files at once
# 'apply' to a 'list' = 'lapply'
# apply 'modified 'read.table''header_fun' function to the list of otu files (otu_path)
# guide: https://nicercode.github.io/guides/repeating-things/
#listapply
ldf<-lapply(otu_path, header_fun)
str(ldf)# this produces a list of dataframes
length(ldf) # length is 4 because there were 4 files imported, and a dataframe made for each file
View(ldf[[2]]) # view second df in list

library(dplyr)
all_otus<-bind_rows(ldf)##combine list of dataframes
View(all_otus)# now they are all in a single file - same result as before

#save
write.csv(all_otus, '/RdataDust/all_otus.csv', row.names=FALSE)

# now this last part (binding the dataframes) will only work if the column headers are exactly the same in each file
# so if there were different numbers of samples between the files...

# melt dataframe so there is a column for each sample name
# basically reduce a subset of columns into 2 columns (1 with the former column name, 1 with the values)
library(reshape2)
View(file2)

file2.melt<-melt(file2, id.vars=c('ID',taxa_list))
View(file2.melt) #compare this to file2

# include syntax to name new columns
file2.melt<-melt(file2, id.vars=c('ID',taxa_list), value.name='count', variable.name='sample')
str(file2.melt)
# to consider another time - 
# how to deal with otus that are same taxomony but different otu?
# are the otu id's 'denovo0' the same otu across files? or these are unique to each file?

## add 'melt' to function
header_melt_fun<-function(x){
  taxa_list<-c('kingdom','phylum','class','order','family','genus','species') #taxonomy column headers
  temp.df<-read.table(otu_path[1], fill=TRUE, stringsAsFactors = FALSE)
  head_list<-as.character(temp.df[1, 1:(length(temp.df)-7)]) #dataframe headers
  out.df<-read.table(x, fill=TRUE, skip=1,col.names=c(head_list,taxa_list)) # read in and assign headers
  out.df$file<-paste(x) # create additional column that designates source file
  out.df.melt<-melt(out.df, id.vars=c('ID','file', taxa_list), value.name='count', variable.name='sample')
  return(out.df.melt)
}

all_otus_melt<-bind_rows(lapply(otu_path, header_melt_fun))
View(all_otus_melt) #everything's here

write.csv(all_otus_melt, 'all_otus_melt.csv', row.names=FALSE)
