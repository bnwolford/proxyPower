#!/usr/bin/env python

#===============================================================================
# Copyright (c) 2019 Brooke Wolford
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

# Python 2.7.6
############################
##### IMPORT MODULES #######
###########################
import subprocess
import argparse
from itertools import islice
import gzip, re, os, math, sys
import copy
import datetime
import numpy as np

###########################
##### PARSE ARGUMENTS ####
###########################
def get_settings():
  parser = argparse.ArgumentParser(description='''Script to perform proxy-case assignment using kinship matrix from KING2, self reported affected relative status of mother, father, sibling, and case/control status from EHR derived phenotypes. The default is an output of the phenotype file with an additional column holding the proxy-case assignment.''')
  parser.add_argument("-k", "--kinship", help="Kinship from KING2 and requires header. Assumes FID1, ID1, FID2, ID2 and additional columns may vary.", type=str)
  parser.add_argument("-ck","--columnKin",help="0-based column number with Kinship value from KING.",default=8,type=int)
  parser.add_argument("-p", "--pheno",help="Tab delimited phenotype file. First column must be an ID specific to the individual sample (e.g. IID).",type=str,required=True)
  parser.add_argument("-d","--header",help="Header",action='store_true')
  parser.add_argument("-o","--output",help="Output file name prefix",type=str,required=True)
  #parser.add_argument("-cm","--columnMother",help="0-based column number for affected mother. Expects 1 if mother is affected and 0 otherwise.", type=int,required=True)
  #parser.add_argument("-cf","--columnFather",help="0-based column number for affected father. Expects 1 if father is affected and 0 otherwise.", type=int,required=True)
  #parser.add_argument("-cs","--columnSibling",help="0-based column number for affected sibling. Expects 1 if sibling is affected and 0 otherwise.", type=int,required=True)
  parser.add_argument("-cp","--columnPhenotype",help="0-based column number for phenotype information. Expects 1 for case, 0 for control, NA for missing [default=12]",type=int,default=12)
  parser.add_argument("-g","--GRS",help="File with ID that matches kinship file and GRS",type=str)
  parser.add_argument("-cg","--columnGRS",help="0-based column number for GRS in -g file",type=int)
                      
  args = parser.parse_args()
  print >> sys.stderr, "%s\n" % args
  return args


############################
######### FUNCTIONS ########
############################


#open file taking zipped or unzipped into account
def openFile(filename):
  if ".gz" in filename:
    command=gzip.open(filename,"rt")
    print >> sys.stderr, "Opening zipped file\n"
  elif ".gz" not in filename:
    command=open(filename,"rt")
    print >> sys.stderr, "Opening unzipped file\n"
  return command


#these numbers are from http://people.virginia.edu/~wc9c/KING/manual.html
def is_first_degree_relative(kinship):
  return (float(kinship) >= 0.177 and float(kinship) <= 0.354)


# read kinship file, assumes FID and IID are equal because no family info
def readKinship(file,col):
  kinDict = {}  # intialize
  openCommand=openFile(file) #handle zipped file
  with openCommand as f:
    next(f)  # skip header
    for line in f:
      line = line.rstrip()
      lineList = line.split("\t")
      IID1=lineList[1]
      IID2=lineList[3]
      KinVal=lineList[col]
      if is_first_degree_relative(KinVal): #record first degree relatives only
        if IID1 not in kinDict.keys():
          kinDict[IID1] = []  # initalize
        kinDict[IID1].append(IID2) #make a list of first degree relative IDs
        if IID2 not in kinDict.keys():
          kinDict[IID2]=[] #initialize
        kinDict[IID2].append(IID1)
        #any ID that is listed in kinship file as ID1  or ID2 will get its own entry in the kinship dictionary
  f.close()
  return kinDict

#read phenotype file with case/control information for sample and affected status of relatives
def readPheno(file,header_bool):
  phenoDict = {}  # initialize
  totalCol=0 #initialize count of columns so we know what is new column to add proxycase assignment
  if header_bool==True:
    count=0
  else:
    count=1
  f = open(file, "r")
  for line in f:
    line = line.rstrip()
    if count==0:
      header=line
      count+=1
    else:
      header=""
      line_list = line.split("\t")
      #replace -9 with NA
      line_list_v2=["NA" if x=="-9" else x for x in line_list]
      phenoDict[line_list[0]] = line_list
      totalCol=len(line_list)
  return phenoDict,totalCol,header

#Perform proxy case assignment using information on case/control status of the sample and the kinship matrix with the rest of samples in the study
def proxy_via_kinship(pd, kd, tc, cp):
  for sample in pd: #for every sample in phenotype file
    if sample in kd: #if sample is in kinship dictionary
      flag=False
      for relative in kd[sample]: #for all listed relatives
        if relative in pd: #if relative is in phenotype file
          if pd[relative][cp]=="1": #if relative is a case
            pd[sample].append("1") #assign positive family history
            flag=True #change flag 
            break
      if flag==False: #if no relatives with disease
        pd[sample].append("0") #assign negative family history
    else: #if sample is not in kinship dictionary
      pd[sample].append("0") #assign negatie family history

##figure out which samples are in the top 5th percentile, return list of those IDS
def percentiles(gd):
  grs_list=[]
  for sample in gd.keys():
    grs_list.append(np.float(gd[sample]))
  p=np.percentile(grs_list,95) #top 5th precentile 
  top_list=[]
  for sample in gd.keys():
    if np.float(gd[sample]) > p:
      top_list.append(sample) #record which samples are in the top 5th percentile
  return(top_list)
  
def match_grs(grs,col,kinDict,out):
  f = open(grs, "r") #open GRS file 
  grsDict={}
  for line in f: #read GRS file into dictionary
    ls = line.rstrip()
    ll=ls.split("\t")
    grsDict[ll[0]]=ll[col]
  top_list=percentiles(grsDict) #get top 5th percentile samples from grs Dict
  o=open(".".join([out,"GRS.txt"]),"w") #open output file
  o2=open(".".join([out,"top5.txt"]),"w") #open output file 2
  o2.write("\t".join(["sample","top5","fracReltop5","numRel"])) #write header
  o2.write("\n")
  for index in kinDict.keys(): #choose the index relative
    sample_list=[]
    score_list=[]
    sample_list.append(index)
    if index in grsDict.keys():
      score_list.append(float(grsDict[index]))
    else:
      score_list.append(np.nan)
    for relative in kinDict[index]: #for relatives in the list corresponding to the index variant 
        sample_list.append(relative)
        if relative in grsDict.keys():
          score_list.append(float(grsDict[relative]))
        else:
          score_list.append(np.nan)
    mean=np.mean(np.array(score_list[1:])) #calculate mean GRS for all relatives of the index
    o.write("\t".join([",".join(str(x) for x in score_list),",".join(sample_list),str(score_list[0]),str(mean)])) #1st value of each list is the index
    o.write("\n")
    
    ##number of relatives in top 5 for index in top 5
    total_relatives=len(sample_list)
    top_relatives=sum(rel in sample_list for rel in top_list)
           
    if index in top_list:
      o2.write("\t".join([index,"1",str(np.float(top_relatives)/np.float(total_relatives)),str(total_relatives)]))
      o2.write("\n")
    else:
      o2.write("\t".join([index,"0",str(np.float(top_relatives)/np.float(total_relatives)),str(total_relatives)]))
      o2.write("\n")
    

  o.close()
  o2.close()
#########################
########## MAIN #########
#########################

def main():
  args = get_settings()
  #always read phenotype file
  phenoDict, totalCol, header = readPheno(args.pheno,args.header)  # read self report file
  print >> sys.stderr, "Finished reading phenotype file %s at %s\n" % (args.pheno, datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
  
  kinDict = readKinship(args.kinship,args.columnKin)  # read kinship file
  print >> sys.stderr, "Finished reading kinship file %s\n" %args.kinship
  
  cp=args.columnPhenotype
  #  cm=args.columnMother
  #  cf=args.columnFather
  #  cs=args.columnSibling

  if args.GRS and args.columnGRS:
    match_grs(args.GRS,args.columnGRS,kinDict,args.output)

    print >> sys.stderr, "Listing GRS per index sample\n"
         
    ############### kinship only ####################

  print >> sys.stderr, "Assigning positive family history based on kinship (-x K) at %s" % datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
  #proxy_via_kinship(phenoDict, kinDict, totalCol,cp) #edits phenodict in place
  print >> sys.stderr, "Finished assigning proxy-case based on kinship at %s\n" % datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')

  f=open(".".join([args.output,"pheno.txt"]),"w")
  if args.header==True:
    header_list=header.split("\t")
    header_list.append("InferredFamHx") #add new column label to header
    f.write("\t".join(header_list))
    f.write("\n")
    for sample in phenoDict:
      f.write("\t".join(phenoDict[sample]))
      f.write("\n")
    f.close()
    print >> sys.stderr, "Finished printing results at %s\n" % datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
  

# call main
if __name__ == "__main__":
  main()

