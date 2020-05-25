###Merve Alp
###2020-05-17
###R studio version 1.1.463
###R version  3.6.2  (2019-12-12)
###version 1
rm(list=ls())

#installing and importing required packages
if   (!require("pacman")) install.packages("pacman")
pacman::p_load("dplyr",
               "tidyr",
               "reshape2",
               "stringr",
               "broom")

library(dplyr)
library(tidyr)
library(reshape2)
library(stringr)
library(broom)

#setting working directory --> where the input csv file is
setwd("~/Desktop/Projects/Other projects/Skyline UGM/") 

#Importing csv file and modifying column names
data              <- read.csv(file= "MSstats_for_Merve_UGM.csv", header=TRUE, sep=",", stringsAsFactors = F, as.is = T, dec = ".")
data$File.Name    <- gsub(pattern = "^[^_]*_", replacement = "",  data$File.Name)###optional ##removes the date (YY.MM.DD_) from the file name
data$File.Name    <- gsub(pattern = "\\..*",   replacement = "",  data$File.Name)##removes .wiff extension
data$File.Name    <- gsub(pattern = "_",       replacement = ".", data$File.Name)###optional  ##Puts dots as seperator instead of underscores
data$Condition    <- gsub(pattern = "_",       replacement = ".", data$Condition)###optional  ##Puts dots as seperator instead of underscores
data$Protein.Name <- gsub(pattern = "^.*\\|",  replacement = "",  data$Protein.Name)###should be modified depending on how protein name look like



dim(data)
names(data)
unique(data$File.Name)

data[c("Area","Background", "Retention.Time", "Start.Time", "End.Time")] <- sapply(data[c("Area",
                                                                                          "Background", 
                                                                                          "Retention.Time", 
                                                                                          "Start.Time", 
                                                                                          "End.Time")], as.numeric)

#data[c("Area","Background")] <- sapply(data[c("Area","Background")], log2)
#data[c("Area","Background")] <- sapply(data[c("Area","Background")], function(x) replace(x, is.infinite(x), NA))


my_zfact <- function(area, background){
        sumsd    <- sd(area, na.rm = TRUE) + sd(background, na.rm = TRUE)
        diffmean <- abs(mean(area, na.rm = TRUE) - mean(background, na.rm = TRUE))
        zfact    <- 1-((3*sumsd) / diffmean)
        return(zfact)
}

test <- data.frame(Area = c(5,5.5,NA), Background = c(1,1.5,2))
my_zfact(test$Area, test$Background)


detect_sym  <- function(start, end, apex ) {
        head <- abs(apex - start)
        tail <- abs(apex - end)
        if(head != 0){
                sym <- (head + tail) / (2* head)  
                return(sym)
        }
        else 
                return(NA)
}

test <- data[15,14:16]
detect_sym(test$Start.Time, test$End.Time, test$Retention.Time)


#Filtering data
data       <- data %>%
        filter(Truncated == "False")   %>% 
        filter(Quantitative == "True") %>%
        filter(Standard.Type != "iRT") %>% 
        select(-c(Truncated, Quantitative, Standard.Type)) 
        

zFact_table <- data %>% filter(Isotope.Label.Type != "light") %>%
        group_by(Protein.Name, Peptide.Modified.Sequence, Fragment.Ion) %>% 
        summarise(z.Factor = my_zfact(Area, Background))  %>%  
        group_by(Protein.Name, Peptide.Modified.Sequence) %>%
        arrange(desc(z.Factor), .by_group = TRUE)

#test <- as.data.frame(data %>%
#                      group_by(Protein.Name, Peptide.Modified.Sequence, Fragment.Ion) )

#test <- data[72,10:12]
peakSym     <- data %>% 
        group_by(Protein.Name, Peptide.Modified.Sequence, Fragment.Ion) %>%
        rowwise()%>%
        summarise(peak.Sym = detect_sym(Start.Time, End.Time, Retention.Time))

peakSym_table <- data.frame(data$Protein.Name, 
                            data$Peptide.Modified.Sequence, 
                            data$Fragment.Ion, 
                            data$Retention.Time, 
                            as.data.frame(peakSym))

#plot(peakSym_table$data.Retention.Time, peakSym_table$peak.Sym)

write.table( zFact_table, file = 'ZfactorResultsRanked.txt',
             quote = F, sep = '\t', dec = '.', row.names = F )


write.table( peakSym_table, file = 'PeakSym.txt',
             quote = F, sep = '\t', dec = '.', row.names = F )
