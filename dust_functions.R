## functions for rdatadust analysis

## split taxa
# a function that splits otu table output with taxonomy listed all in one column to seprate columns for each grouping
# provide dataframe and name of column containing the taxonomy list
# 

split_taxa<-function(df, col='taxonomy'){
  df$kingdom<-str_match(df[,col], "k__(.*?);")[,2]
  df$phylum<-str_match(df[,col], "p__(.*?);")[,2]
  df$class<-str_match(df[,col], "c__(.*?);")[,2]
  df$order<-str_match(df[,col], "o__(.*?);")[,2]
  df$family<-str_match(df[,col], "f__(.*?);")[,2]
  df$genus<-str_match(df[,col], "g__(.*?);")[,2]
  df$species<-str_match(df[,col], "s__(.*?)")[,2]
  return(df)
}

# caluclate standard error (unaffected by na's)
se <- function(x) sqrt(var(x)/length(x))
