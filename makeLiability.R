#!/usr/bin/Rscript

#Revised from make_liabs.R script by Jimmy Liu (https://github.com/jimmyzliu/pmliability)

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
    make_option(c("-t", "--heritability"), type="numeric", help="Estimated eritability of the trait"),
    make_option(c("-v","--prevalence"), type="character", help="File with prevalence by age and sex"),
    make_option(c("-p","--ped"),type="character",help=".ped file"),
    make_option(c("-o","--output"),type="character",help="Output file name"),
    make_option(c("-f","--phenoFile"),type="character",help="Phenotype file. If present, will add liabilities as a new column matching on IID")
	        )

parser <- OptionParser(
    usage="%prog -t heritability -v prevalence -p ped_file -o output", add_help_option=TRUE,
        option_list=optionList
	    )

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

########################### Libraries ###################################################

library(ggplot2)
library(tmvtnorm) #https://cran.r-project.org/web/packages/tmvtnorm/tmvtnorm.pdf
library(kinship2)  #https://cran.r-project.org/web/packages/kinship2/kinship2.pdf
library(data.table)
library(dplyr)
library(GENESIS) #bioconductor package to turn king2 output into matrix
library(Matrix)

########################### Functions ################################

#if age from .ped is an outlier and not present in prevalence file, add a year until a prevalence is available or until max age is reached 
#if age is missing use median proband age 
assign_age<-function(fam_info,i,prev,med_age){
	age <- fam_info$AGE[i]
	if(is.na(age)){ age <- med_age }
	repeat{
		if (age>max(prev$AGE)){
		   break
		} else if (nrow(prev[prev$AGE==age])==1){
		   break
		}
		age <- age + 1
	   }
	 return(age)
}


read_matrix<-function(mat){
    m<-king2mat(mat,type="kinship") #convert KING2 to m
    #sm<-as(m,'sparseMatrix') #make dgCMatrix
    sm<-Matrix(m,sparse=TRUE) #make dsCmatrix
                                        #maybe can make faster with Rcpp http <- //gallery.rcpp.org/articles/sparse-matrix-coercion/
    
    }


############################ Main Commands ##################################################

h2 <- opt$heritability
pedigree <-opt$ped
prev<- opt$prevalence
out <- opt$output
phenoFile <-opt$phenoFile

#check for required arguments
if (is.null(h2) | is.null(pedigree) | is.null(prev) | is.null(out)) {
   stop("Missing  argument!\n")
}

prev <- fread(prev,header=T)
ped <- fread(pedigree,header=T)
fams <- unique(ped$FID)

# 1 == female, 0 == male
ped$FATHER[which(is.na(ped$FATHER))] <- ""
ped$MOTHER[which(is.na(ped$MOTHER))] <- ""

# change sex to: 1 == male, 2 == female, 3 == unknown
ped$SEX[which(is.na(ped$SEX))] <- "unknown" #could also be -9?
ped$SEX[which(ped$SEX == "2")] <- "female"
ped$SEX[which(ped$SEX == "1")] <- "male"

#must have both mom and dad in ped
#big_ped <- pedigree(id = ped$IID,
 #                   dadid = ped$FATHER,
  #                  momid = ped$MOTHER,
   #                 sex = ped$SEX,
    #                famid = ped$FID,
     #               affected = ped$PHENO)

#calculate median age of probands
med_age<-median(ped[ped[,lapply(.SD,function(x) x %like% "PROBAND")]$IID,]$AGE)

liab_pheno <- data.frame()
for(fam in fams){
    fam_info <- ped[ped$FID == fam,]
    fam_ped <- pedigree(id = fam_info$IID, #create pedigree structure for each family
                        dadid = fam_info$FATHER,
                        momid = fam_info$MOTHER,
                        sex = fam_info$SEX,
                        famid = fam_info$FID,
                        affected = fam_info$PHENO)
    #fam_ped <- big_ped[fam]
    sigma <- 2*kinship(fam_ped) #matrix of kinship coefficients from family pedigree x 2 = coefficient of relationship/genetic relatedness
    sigma <- h2 * sigma #sigma is dsCMatrix
    diag(sigma) <- rep(1,nrow(sigma)) #diaganol is 1 
    ids <- fam_ped$id

    #if(NA %in% fam_info$PHENO){
     #next
     #}
    prevs <- rep(0,length(ids))

    #upper
    l_upper <- rep(0,nrow(fam_info))
    l_lower <- rep(0,nrow(fam_info))
    for(i in 1:nrow(fam_info)){
    	age <- assign_age(fam_info,i,prev,med_age) #deal with corner cases regarding missing or outlier age in pedigree
	sex <- fam_info$SEX[i]
        pheno <- fam_info$PHENO[i]
        # if missing phenotypes in family, exclude. fix later
        if(age < min(prev$AGE)){ age <- min(prev$AGE)}
        if(age > max(prev$AGE)){ age <- max(prev$AGE)}
        k <- sum(prev[prev$AGE == age,c(2,3)])/2 #average male and female prevalence
        if(sex == "male"){ k <- as.numeric(prev[prev$AGE == age,4])} #male smooth prevalence
        if(sex == "female"){ k <- as.numeric(prev[prev$AGE == age,5])} #female smooth prevalence

	#k cannot be negative
        t <- qnorm(1-k)
        if(is.na(pheno)){
            l_upper[i] <- Inf
            l_lower[i] <- -Inf
        }
        else{
            if(pheno == 1){ #case
                l_upper[i] <- Inf
                l_lower[i] <- t
                if(k == 0){
                    all_prev <- unlist(c(prev[,2],prev[,3]))
                    min_k <- min(all_prev[all_prev!=0])
                    t <- qnorm(1-min_k)
                    l_lower[i] <- t
                }
            }
            if(pheno == 0){ #control
                l_upper[i] <- t
                l_lower[i] <- -Inf
            }
        }
    }
    liab <- rtmvnorm(1000,
                     mean = rep(0,length(ids)),
                     sigma = sigma,
                     lower = l_lower,
                     upper = l_upper,
                     algorithm = "gibbs",
                     burn.in.samples = 100,
                     thinning = 5) #decide optimal thinnign with autocorrelation plots
    liab <- colSums(liab)/1000 
    liab_pheno <- rbind(liab_pheno,data.frame(ids,liab))
    if(nrow(liab_pheno) %% 500 == 0){ print(nrow(liab_pheno))} #print line number every 500 lines to show progress
}


#write table of results
file<-paste0(out,".liab")
write.table(liab_pheno,file,col.names=F,row.names=F,quote=F)

#make plot
pdf_name<-paste0(out,".liab.pdf")
pdf(file=pdf_name,height=3,width=3)
title=paste(sep=" ",out,"Liabilities")
#ggplot(liab_pheno,aes(x=liab)) + theme_bw() + labs(x="Posterior mean liability") + geom_density()
ggplot(liab_pheno,aes(x=liab)) + theme_bw() + labs(x="Posterior mean liability",y="Frequency",title=title) + geom_histogram(binwidth=0.01)
dev.off()

#add liability of probands back to phenotype file if option is given
if (!is.null(phenoFile)) {
    p<-fread(phenoFile,header=T) #read in phenotype file
    probands<-liab_pheno[grepl('PROBAND',liab_pheno$ids),] #subset to probands only
    probands$new_id<-gsub('_PROBAND',"",probands$ids)
    tmp<-left_join(p,probands,by=c("IID"="new_id")) #join the calculated liabilities to the phenotype file
    new_p<-tmp[,!names(tmp) %in% c("ids")]

    ##write file
    file<-paste0(out,".liab.phenoFile")
    write.table(new_p,file,col.names=T, row.names=F, quote=F,sep="\t")
    
    ## plot just the proband liabilities 
    pdf_name<-paste0(file,".pdf")
    pdf(file=pdf_name,height=3,width=3)
    title=paste(sep=" ",out,"Proband Liabilities"))
    ggplot(new_p,aes(x=liab)) + theme_bw() + labs(x="Posterior mean liability",y="Frequency",title=title) + geom_histogram(binwidth=0.01)
    dev.off()
    
}
