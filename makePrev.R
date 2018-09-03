#!/usr/bin/Rscript

#could use rprev eventually (https://cran.r-project.org/web/packages/rprev/vignettes/user_guide.html)

#===============================================================================
# Copyright (c) 2018 Brooke Wolford
# Lab of Dr. Cristen Willer and Dr. Mike Boehnke
# University of Michigan

#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.
#==============================================================================


################### Read Command Line Parameters #############################
suppressPackageStartupMessages(library(optparse))
optionList <- list(
    make_option(c("-f", "--file"), type="character", help="Phenotype file with header"),
    make_option(c("-b","--birthYear"),type="character",help="Column name of birthYear",default="BirthYear"),
    make_option(c("-s","--sex"),type="character",help="Column name of sex",default="Sex"),
    make_option(c("-t","--trait"),type="character",help="Column name of trait"),
    make_option(c("-o","--output"),type="character",help="Output file name"),
    make_option(c("-u","--subset"),action="store_true",help="This flag subsets phenotype file to persons born from 1940-1960",default="False")
)

parser <- OptionParser(
    usage="%prog -f <file>",
    option_list=optionList
)

# a hack to fix a bug in optparse that won't let you use positional args
# if you also have non-boolean optional args:
getOptionStrings <- function(parserObj) {
        optionStrings <- character()
        for(item in parserObj@options) {
            optionStrings <- append(optionStrings, c(item@short_flag, item@long_flag))
        }
        optionStrings
}

optStrings <- getOptionStrings(parser)
arguments <- parse_args(parser, positional_arguments=TRUE)

#exit program if parameters are not provided
if (is.na(arguments$options$file)){
    quit()
}

########################### Load functions and libraries ###################################################

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

################################  DEFINE FUNCTIONS ###################################


###################################### Main commands $###############################

#### file names #####
file <- arguments$options$file
cb<- arguments$options$birthYear
cs<-arguments$options$sex
ct<-arguments$options$trait
out <-arguments$options$output
subset<-arguments$options$subset

#read in file
if (grepl('.gz',file)) {
    df<-fread(paste(sep=" ","zcat",file),header=TRUE)
} else {
    df <- fread(file,header=TRUE)
}

#subset to everyone born 40s-60s, change this to be an option later
if (subset) {
   age_filter<-paste0(cb,">=1940 & ", cb, "<=1960") #work around for filter
   sub<-as_tibble(df) %>% filter_(age_filter)
} else {
  sub<-df
}
sub$Age<-2018-sub$BirthYear #assumes current year is 2018

total<-sub %>% group_by_("Age",cs) %>% tally() %>% spread_(cs,"n")
trait_def=paste0(ct,"==2") #work around for filter 
trait<-sub %>% group_by_("Age",cs,ct) %>% tally() %>% filter_(trait_def) %>% spread_(cs,"n") %>% select(-one_of(ct))
names(trait)<-c("Age","Male","Female") #this assumes Male is 1 and Female is 2 in phenotype file
names(total)<-c("Age","Male","Female")

#calculate prevalence from counts of total and trait incidence per age/sex category
prev<-full_join(total,trait,by="Age")
prev<-prev %>% mutate(Male=Male.y/Male.x) %>% mutate(Female=Female.y/Female.x) %>%  mutate_if(is.numeric, funs(replace(., is.na(.), 0))) %>% select(Age,Male,Female)

#plot
plot_df<-as_tibble(prev) %>% gather(key=Sex,value=Prev,Male,Female) #reformat wide to long for plotting
plot_name<-paste(sep=".",out,"pdf")
title=paste0("Prevalence of ",ct)
pdf(file=plot_name,height=3,width=5)
ggplot(plot_df,aes(x=Age,y=Prev,color=Sex,shape=Sex)) + geom_point()  + theme_bw() + scale_color_manual(values=c("red","blue")) + labs(title=title,xlab="Age",ylab="Prevalence(K)") + geom_smooth(method="loess") + scale_y_continuous(limits = c(0, NA)) + scale_x_continuous(limits=c(NA,100))
dev.off()

#save Loess smoothed values
prev$MaleSmooth<-predict(loess(Male~Age,prev),prev$Age)
prev$FemaleSmooth<-predict(loess(Female~Age,prev),prev$Age)

#restrain predictions to (0,1) boundaries
prev[prev$MaleSmooth<0,]<-0
prev[prev$MaleSmooth>1,]<-1
prev[prev$FemaleSmooth<0,]<-0
prev[prev$FemaleSmooth>1,]<-1

#rename to match makeLiability.R
names(prev)<-c("AGE", "K_MALE", "K_FEMALE", "K_MALE_SMOOTH","K_FEMALE_SMOOTH")

#write out file
dat_name<-paste(sep=".",out,"tab")
write.table(prev,append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE, file=dat_name)