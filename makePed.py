#!/usr/bin/env python

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
import random

###########################
##### PARSE ARGUMENTS ####
###########################
def get_settings():
  parser=argparse.ArgumentParser(
  description='''Script to create .ped file from phenotype file. Creates dummy F and M entries from self report affected family member''')
  parser.add_argument("-p","--pheno",help="Phenotype file", type=str, required=True)
  parser.add_argument("-cm","--columnMother",help="0-based column number for affected mother. Expects 1 if mother is affected and 0 otherwise.", type=int)
  parser.add_argument("-cf","--columnFather",help="0-based column number for affected father. Expects 1 if father is affected and 0 otherwise.", type=int)
  parser.add_argument("-cs","--columnSibling",help="0-based column number for affected sibling. Expects 1 if sibling is affected and 0 otherwise.", type=int)
  parser.add_argument("-cp","--columnPhenotype",help="0-based column number for phenotype information. Expects 1 for case, 0 for control, NA for missing [default=12]",type=int,default=12)
  parser.add_argument("-cr","--columnRelative",help="0-based column number of any affected first degree relative. Expects 1 if a 1st degree relative is affected.", type=int)
  parser.add_argument("-o","--outputFile",help="Prefix for output ped file.",type=str,required=True)
  parser.set_defaults(remove=False)
  args=parser.parse_args()
  return args
  sys.stdout.write(args)
  sys.stdout.write(sys.version)

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

#make pedigree if you have info on affected mother, father, sibling
def makePed_fam(phenoDict,cm,cf,cs,cp):

  return 0

#make pedigree if you have info on affected first degree relative in general
# in output 1 is case and 0 is control and NA is missing
# 2 is female and 1 is male
def makePed(phenoDict,cp,cr,out):
  for id,value in phenoDict.iteritems():
    #figure out phenotype of proband
    proband_id="_".join([value[0],"PROBAND"])

    #randomly choose sex of affected parent
    parent=random.randint(1,2)
    mom_id="_".join([value[0],"MOTHER"])
    dad_id="_".join([value[0],"FATHER"])
    parent_entry="\t".join([dad_id,mom_id])
    if value[11]=="2": #case
      out.write("\t".join([value[0],proband_id,parent_entry,str(value[4]),str(2018-int(value[5])),str(1)])) #print proband line
    elif value[11]=="1": #control
      out.write("\t".join([value[0],proband_id,parent_entry,str(value[4]),str(2018-int(value[5])),str(0)])) #print proband line
    else: #NA
      out.write("\t".join([value[0],proband_id,parent_entry,str(value[4]),str(2018-int(value[5])),"NA"])) #print proband line
    out.write("\n")
    
    #create dummy mom and dad
    parent_age=str(2018-int(value[5])+20) #estimating parents to be 20 years older
    if value[12]=="2": #yes affected relative
      if parent==2: #make mom affected
        out.write("\t".join([value[0],mom_id,"NA","NA","2",parent_age,str(1)]))
        out.write("\n")
        out.write("\t".join([value[0],dad_id,"NA","NA","1",parent_age,str(0)]))
      elif parent==1: #make dad affected
        out.write("\t".join([value[0],mom_id,"NA","NA","2",parent_age,str(0)]))
        out.write("\n")
        out.write("\t".join([value[0],dad_id,"NA","NA","1",parent_age,str(1)]))
    elif value[12]=="1": #no affected relative
      out.write("\t".join([value[0],mom_id,"NA","NA","2",parent_age,str(0)]))
      out.write("\n")
      out.write("\t".join([value[0],dad_id,"NA","NA","1",parent_age,str(0)]))
    else: #missing or unknown
      out.write("\t".join([value[0],mom_id,"NA","NA","2",parent_age,"NA"]))
      out.write("\n")
      out.write("\t".join([value[0],dad_id,"NA","NA","1",parent_age,"NA"]))
    out.write("\n")
      
  return 0
  
#########################
########## MAIN #########
#########################

def main():  
    args = get_settings()

    phenoDict,header=readPheno(args.pheno) #read phenotype file with proxy-case definition

    out=".".join([args.outputFile,"ped"])
    o=open(out,"w")
    o.write("\t".join(["FID","IID","FATHER","MOTHER", "SEX ","AGE", "PHENO\n"])) #should we be working in the birthYear space?
    
    random.seed(12345) #set seed 
    makePed(phenoDict,args.columnPhenotype,args.columnRelative,o)

    o.close()
    
#call main
if __name__ == "__main__":
  main()

