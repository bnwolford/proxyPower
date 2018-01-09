#!/usr/bin/env python

#===============================================================================
# Copyright (c) 2017 Brooke Wolford
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

#Python 2.7.6
############################
##### IMPORT MODULES #######
###########################
import subprocess
import argparse
#import numpy
#import rpy2
from itertools import islice
import gzip, re, os, math, sys
import copy

###########################
##### PARSE ARGUMENTS ####
###########################
def get_settings():
  parser=argparse.ArgumentParser(
  description='''Script to convert proxy-case assignment phenotype file to a phenotype file ready for analysis. Model definitions are 1=standard GWAS, 2=GWAS with proxy-cases removed from controls (i.e. cleaner controls), 3=GWAX with proxy-cases as cases, 4=Cases vs proxy-cases vs controls (i.e. appropriately modelling proxy-cases as intermediate), 5=Cases + proxy-cases vs controls. Identical to --pheno file except --column has been converted to values that correspond to the --model.'''
)
  parser.add_argument("-k","--kinship",help="Kinship from KING2 and requires header", type=str,required=True)
  parser.add_argument("-p","--pheno",help="Phenotype file of any format. For Models 2-5 the F column must be 0, 0.25, 0.5, 1 or NA.\nDefault expects format IID FID PATID MATID Sex BirthYear batch PC1 PC2 PC3 PC4 F. Requires header line. This file is read into memory.", type=str, required=True)
  parser.add_argument("-c","--column",help="0-based column number for F column with 0, 0.25, 0.5, 1 or NA [default=12]",type=int,default=11)
  parser.add_argument("-o","--output",help="Full path for name and location of output file. Output file is ready for BOLT-LMM with --phenoCol=F.", type=str, required=True)
  parser.add_argument("-m","--model",help="Type of model and way to consider proxy-cases\n[1=standard GWAS, 2=GWAS with cleaner controls, 3=GWAX, 4=Cases vs proxy-cases vs controls, 5=Cases + proxy-cases vs controls]\n", type=int, required=True)
  parser.add_argument("-r","--remove2dr",help="Use flag to print second degree relatives (F=0.25) as NA thereby removing those samples from analysis [default=FALSE]",action="store_true",dest='remove')
  parser.set_defaults(remove=False)
  args=parser.parse_args()
  return args

############################
######### FUNCTIONS ########
############################

#read in phenotype file
def readPheno(file):
  phenoDict={} #initialize
  count=0
  f=open(file,"r")
  for line in f:
    line = line.rstrip()
    if count==0:
      header=line #save header
      count+=1
    else:
      line_list=line.split("\t")
      phenoDict[line_list[0]]=line_list #ID is key and full line is value
      count+=1
  return phenoDict, header


#creates phenotype file for BOLT-LMM that uses the proxy-cases as the specified model dictates
def updateF(phenoDict,col,model,output,header):

  #delete output file if it exists because we are appending
  try:
    os.remove(output)
  except OSError:
    pass
  f1 = open(output, 'a')

  print >> f1, header



  for line in phenoDict.keys(): #loop through every line in phenotype file
    
    zc=col #0-based column for F
    f=str(phenoDict[line][zc]) #save F value as string

    if model == 1: #standard GWAS
      print >> sys.stderr, "Model 1 is Standard GWAS and a phenotype file for this analysis is created by proxyCaseAssign1dr.py or proxyCaseAssignAffRel.py based on the phenotype files provided there.\n"
      break

    elif model == 2: #GWAS with cleaner controls
      if f=="0.25":
        phenoDict[line][zc]="NA" #convert 2nd degree proxy-case to NA
      elif f=="0.5":
        phenoDict[line][zc]="NA" #convert 1st degree proxy-case to NA

    elif model == 3: #GWAX using proxy-cases as cases
      if f=="0.25":
        phenoDict[line][zc]="NA" #convert 2nd degree proxy-case to NA
      elif f=="1": #convert cases to NA
         phenoDict[line][zc]="NA" #convert cases to NA
      elif f=="0.5":
        phenoDict[line][zc]=1 #convert 1st degree proxy-case to case 

    elif model == 4: #model cases, proxy-cases, and controls separately as semi-continuous
      if f=="0.25":
        phenoDict[line][zc]="NA" #convert 2nd degree proxy-case to NA

    elif model == 5: #GWAS grouping proxy-cases and cases as cases 
      if f=="0.25":
        phenoDict[line][zc]="NA" #convert 2nd degree proxy-case to NA
      elif f=="0.5":
        phenoDict[line][zc]=1 #convert 1st degree proxy-case to case

    else:
      print >> sys.stderr, "Model variable is not expected. Please enter 1, 2, 3 or 4.\n"

    print >> f1, "\t".join(str(x) for x in phenoDict[line]) #print to output file
  
  f1.close() #close file  

def updateF_2dr(phenoDict,col,model,output,header):
  print >> sys.stderr, "Functionality not yet available.\n"

#########################
########## MAIN #########
#########################

def main():  
    args = get_settings()

    phenoDict,header=readPheno(args.pheno) #read phenotype file with proxy-case definition
    
    if args.remove:
      updateF(phenoDict,args.column,args.model,args.output,header) #update F to match model
    elif args.remove == False: 
      #updateF_2dr(phenoDict,args.column,args.model,args.output,header) #update F to match model, keep second degree relatives
      print >> sys.stderr, "Functionality not yet available.\n"

#call main
if __name__ == "__main__":
  main()

