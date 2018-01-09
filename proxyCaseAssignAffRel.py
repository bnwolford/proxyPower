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

# Python 2.7.6
############################
##### IMPORT MODULES #######
###########################
import subprocess
import argparse
from itertools import islice
import gzip, re, os, math, sys
import copy


###########################
##### PARSE ARGUMENTS ####
###########################
def get_settings():
  parser = argparse.ArgumentParser(description='''Script to perform proxy-case assignment using kinship matrix from KING2, self reported affected relative status of mother, father, sibling, and case/control status from EHR derived phenotypes. The default is an output of the phenotype file with an additional column holding the proxy-case assignment.''')
  parser.add_argument("-k", "--kinship", help="Kinship from KING2 and requires header", type=str)
  parser.add_argument("-p", "--pheno",help="Tab delimited phenotype file. First column must be an ID specific to the individual sample (e.g. IID). Header expected",type=str,required=True)
  parser.add_argument("-cm","--columnMother",help="0-based column number for affected mother. Expects 1 if mother is affected and 0 otherwise.", type=int,required=True)
  parser.add_argument("-cf","--columnFather",help="0-based column number for affected father. Expects 1 if father is affected and 0 otherwise.", type=int,required=True)
  parser.add_argument("-cs","--columnSibling",help="0-based column number for affected sibling. Expects 1 if sibling is affected and 0 otherwise.", type=int,required=True)
  parser.add_argument("-cp","--columnPhenotype",help="0-based column number for phenotype information. Expects 1 for case, 0 for control, NA for missing [default=12]",type=int,default=12)
  parser.add_argument("-o", "--output",help="Type of output file (BOLT-LMM=B, PLINK=P); default is additional column at end of phenotype file", type=str)
  parser.add_argument("-cc", "--conservControl", help="Requires conservative control (unknown or missing on self report questions are not allowed to be controls), less stringent control is default",action="store_true")
  parser.add_argument("-x", "--proxy",help="Type of logic used to identify proxy-cases [all=A, kinship only=K, self report only=SR, self report minus kinship=SMK, self report plus kinship=SPK]",type=str,required=True)
  parser.add_argument("-n","--number",help="Name of file in which to print number of cases/proxy-cases every sample is related to. If file is not provided then this functionality will not happen.",type=str)
  parser.add_argument("-m","--model1",help="Name of file in which to print model 1 (standard GWAS) phenotype file. This file will be similar to --pheno input file, but header for phenotype column will be F and values will be 1 for case, 0 for control, NA for missing for consistency with proxyModel.py phenotype files\n",type=str,required=True)
  args = parser.parse_args()
  return args


############################
######### FUNCTIONS ########
############################

#these numbers are from http://people.virginia.edu/~wc9c/KING/manual.html
def is_first_degree_relative(kinship):
  return (float(kinship) >= 0.177 and float(kinship) <= 0.354)


# read kinship file, assumes FID and IID are equal because no family info
def readKinship(file):
  kinDict = {}  # intialize
  f = open(file, "r")
  next(f)  # skip header
  for line in f:
    line = line.rstrip()
    (FID1, IID1, FID2, IID2, NSNP, HETHET, IBS0, KinVal) = line.split("\t")
    if IID1 not in kinDict.keys():
      kinDict[IID1] = {}  # initalize
    kinDict[IID1][IID2] = KinVal
  return kinDict

#read phenotype file with case/control information for sample and affected status of relatives
def readPheno(file):
  phenoDict = {}  # initialize
  totalCol=0 #initialize count of columns so we know what is new column to add proxycase assignment
  count=0
  f = open(file, "r")
  for line in f:
    line = line.rstrip()
    if count==0:
      header=line
      count+=1
    else:
      line_list = line.split("\t")
      phenoDict[line_list[0]] = line_list
      totalCol=len(line_list)
  return phenoDict,totalCol,header

#Perform proxy case assignment using information on case/control status of the sample and the kinship matrix with the rest of samples in the study
def proxy_via_kinship(pd, kd, cc, tc, cp):
  pd_kinship = copy.deepcopy(pd)

  if cc:
    print >> sys.stderr, 'Conservative control functionality is not available with kinship only option\n'

  for sample in pd:
    if pd[sample][cp] == '1':  # if case
      pd_kinship[sample].append("1")  # set as case
    elif pd[sample][cp] == '0':  # if control
      pd_kinship[sample].append("0")  # set as control
    else: #if missing or NA
      pd_kinship[sample].append("NA") #missing
      #sys.stderr.write("Sample %s is neither case (1) nor control (0).\n" % sample)

  for ID1, v in kd.items():  # for every ID1 in kinship dictionary
    for ID2, kinship in v.items():  # for every ID2 in kinship dictionary
      if ID1 in pd.keys() and ID2 in pd.keys():

        if is_first_degree_relative(kinship):  # if first degree relative
          if pd[ID1][cp] == '1' and pd[ID2][cp] == '0':  # if case and unaffected
            if pd_kinship[ID2][tc] == 'NA' or float(pd_kinship[ID2][tc]) < 0.5: #if control or missing
              pd_kinship[ID2][tc] = "0.5"  # assign as proxy-case
          elif pd[ID1][cp] == '0' and pd[ID2][cp] == '1':  # if unaffected and case
            if pd_kinship[ID1][tc] == 'NA' or float(pd_kinship[ID1][tc]) < 0.5: #if control or missing
              pd_kinship[ID1][tc] = "0.5"  # assign as proxy-case
        #else:  
          #print >> sys.stderr, "%s and %s are not first degree relatives based on kinship\n"  % (ID1,ID2)

  return pd_kinship

#Perform proxy case assignment using information on case/control status of the sample and the self reported status of mother, father, sibling
def proxy_via_selfreport(pd, cc, cp, cm, cf, cs):
  pd_sr = copy.deepcopy(pd)

  for sample in pd:  # for every sample id in phenotype file
    if pd[sample][cp] == '1':  # if case
      pd_sr[sample].append("1")  # set as case

    elif pd[sample][cp] == '0':  # if unaffected

      #To Do: add functionality so we can define proxy-case based on granularity of WHO is affected (parent vs sibling)
      if pd[sample][cm] == '1' or pd[sample][cf] == "1" or pd[sample][cs] == "1" :  # if have an affected first degree relative
        pd_sr[sample].append("0.5")  # set as proxy-case

      elif pd[sample][cm] == '0' and pd[sample][cf] == "0" and pd[sample][cs] == "0":  # if dont have an affected first degree relative
        pd_sr[sample].append("0") #set as control

      else: #if missing or NA affected first degree relative status
        if cc: #if conservative control option
          pd_sr[sample].append("NA")
        else:
          pd_sr[sample].append("0") #set as control if conservative control option not invoked
    
    else: #if missing or NA for case
      pd_sr[sample].append("NA") #missing
      #sys.stderr.write("Sample %s is neither case (1) nor control (0)\n." % sample)
      
  return pd_sr

#Refine proxy case assignment by self report by considering kinship matrix (e.g. if a proxy-case's affected relative is a case in the study the proxy-case will become NA)
def proxy_via_selfreport_minus_kinship(pd, kd, tc):
  pd_smk = copy.deepcopy(pd)

  for ID1, v in kd.items():  # for every ID1 in kinship dictionary
    for ID2, kinship in v.items():  # for every ID2 in kinship dictionary
      if is_first_degree_relative(kinship):  # if first degree relative
        if ID1 in pd.keys() and ID2 in pd.keys():  # if IDs from kinship matrix are in phenotype file
          # reassign to NA if proxy-case or 2nd degree proxy-case has case in cohort
          if (pd[ID1][tc] == "0.5" and pd[ID2][tc] == "1"):
            pd_smk[ID1][tc] = "NA"
          elif (pd[ID1][tc] == "1" and pd[ID2][tc] == "0.5"):
            pd_smk[ID2][tc] = "NA"
          # reassign to NA if control has case in cohort
          if (pd[ID1][tc] == "0" and pd[ID2][tc] == "1"):
            pd_smk[ID1][tc] = "NA"
          elif (pd[ID1][tc] == "1" and pd[ID2][tc] == "0"):
            pd_smk[ID2][tc] = "NA"
      #else:  
        #print >> sys.stderr, "%s and %s are not first degree relatives based on kinship, leaving proxy-case assignment as is from self report\n" % (ID1,ID2)

  return pd_smk

#Refine proxy case assignemnt by self report by considering kinship matrix (e.g. if a control does not report an affected first degree relative but we identify one in the study using kinship matrix, the control will become a proxy-case)
def proxy_via_selfreport_plus_kinship(pd, kd, tc):
  pd_spk = copy.deepcopy(pd)

  for ID1, v in kd.items():  # for every ID1 in kinship dictionary
    for ID2, kinship in v.items():  # for every ID2 in kinship dictionary

      if is_first_degree_relative(kinship):  # if first degree relative
        if ID1 in pd.keys() and ID2 in pd.keys():  # if IDs from kinship matrix are in phenotype file
          # reassign to 0.5 if control is 1dr to a case or NA from self report (consv control) is related to a case
          if (pd[ID1][tc] == "0" and pd[ID2][tc] == "1") or (pd[ID1][tc] == "NA" and pd[ID2][tc] == "1"):
            if pd_spk[ID1][tc] == 'NA' or float(pd_spk[ID1][tc]) < 0.5:
              pd_spk[ID1][tc] = "0.5"
          elif (pd[ID1][tc] == "1" and pd[ID2][tc] == "0") or (pd[ID1][tc] == "1" and pd[ID2][tc] == "NA"):
            if pd_spk[ID2][tc] == 'NA' or float(pd_spk[ID2][tc]) < 0.5:
              pd_spk[ID2][tc] = "0.5"

      #else:  
          #print >> sys.stderr, "%s and %s are not first degree relatives based on kinship, leaving proxy-case assignment as is from self report\n" % (ID1,ID2)

  return pd_spk

#print output formatted for BOLT-LMM, requires first 10 columns of phenoFile are "FID", "IID", "PATID", "MATID", "Sex", "BirthYear", "batch", "PC1", "PC2", "PC3", "PC4" and last column is F which holds proxy-case assignment
def BOLT_print(pd, tc):
  # assumes BOLT-LMM sees -9 and -NA as missing data in --phenoFile (--phenoCol will be F)
  header = ["FID", "IID", "PATID", "MATID", "Sex", "BirthYear", "batch", "PC1", "PC2", "PC3", "PC4", "F"]  # list of header strings
  print "\t".join(header)  # print header
  for sample in pd:
    data = []
    for x in pd[sample][0:11]: #assumes first 10 columns are as seen in header
      data.append(x)
    data.append(pd[sample][tc]) #F column which holds proxy case assignment
    print "\t".join(str(x) for x in data)

#if -n is provided an output file name, count the number of cases and proxy-cases a given sample is related to within the study
def count_relatives(file,kd,pd,tc):

  f1 = open(file, 'w+')

  for sample in pd: #for every sample in phenotype file

    list=search_nested_dict(kd,sample) #find all relatives

    #initialize counts
    proxy_case_count=0
    case_count=0
    relative_count=len(list)

    #identify if relative is proxy case or case
    for relative in list:
      if relative in pd.keys(): #if relative is in phenotype file
        if pd[relative][tc]=="0.5":
          proxy_case_count+=1
        elif pd[relative][tc]=="1":
          case_count+=1

    data_to_print=[sample,proxy_case_count,case_count,relative_count]
    print >> f1, "\t".join(str(x) for x in data_to_print)

  f1.close()

def search_nested_dict(dict,key):
  found=[]
  for masterKey,masterValue in dict.iteritems():
    if key in masterKey:
      for ID, v in masterValue.items():
        found.append(ID)
    for childKey,childValue in masterValue.iteritems():
      if key in childKey:
          found.append(masterKey)
  return found


#print model 1 (standard gwas) phenotype file based on --pheno but consistent with proxyModel.py output
def model1_print(header,cp,pd,name):
  #delete output file if it exists because we are appending
  try:
    os.remove(name)
  except OSError:
    pass
  f1 = open(name, 'a')

  header_list=header.split("\t")
  header_list[cp]="F" #replace header label with F

  print >> f1, "\t".join(header_list)

  #expects phenotype column to have 1 for case, 0 for control, NA for missing so print as is
  for sample in pd:
    print >> f1, "\t".join(pd[sample])

  f1.close() #close file
  return

#########################
########## MAIN #########
#########################

def main():
  args = get_settings()

  try:
    args.number
  except NameError:
    args.number = None

  #always read phenotype file
  phenoDict, totalCol, header = readPheno(args.pheno)  # read self report file

  print >> sys.stderr, "Finished reading phenotype file %s\n" % args.pheno
  
  if (args.number is not None) or args.proxy=="SMK" or args.proxy=="SPK" or args.proxy=="A" or args.proxy=="K": #only read kinship file into memory if you have to
    kinDict = readKinship(args.kinship)  # read kinship file
    print >> sys.stderr, "Finished reading kinship file %s\n" %args.kinship

  # print kinDict
  # print phenoDict

  cp=args.columnPhenotype
  cm=args.columnMother
  cf=args.columnFather
  cs=args.columnSibling

  # create model 1 (standard gwas) phenotype file
  model1_print(header,cp,phenoDict,args.model1)
  print >> sys.stderr, "Finished printing model 1 phenotype file %s\n" % args.model1

  ########### Self report minus kinship ################
  if args.proxy == "SMK":  # self report minus kinship
    print >> sys.stderr, "Assigning proxy-cases based on self report (-x SMK)"

    phenoDict_SR = proxy_via_selfreport(phenoDict, args.conservControl, cp, cm, cf, cs)  # self report
    print >> sys.stderr, "Finished assigning proxy-cases based on self report\n"
    phenoDict_SRMK = proxy_via_selfreport_minus_kinship(phenoDict_SR, kinDict, totalCol)  # self report minus kinship
    print >> sys.stderr, "Finished refining proxy-case assignment based on kinship\n"

    if args.number is not None:
      count_relatives(args.number,kinDict,phenoDict_SRMK,totalCol)
      print >> sys.stderr, "Finished counting relatives of each sample who are proxy-cases or cases\n"

    if args.output == "B":
      BOLT_print(phenoDict_SRMK,totalCol)
    elif args.output == "P":
      print >> sys.stderr, "This functionality not available yet\n"
    else:
      header_list=header.split("\t")
      header_list.append("F") #add new column label to header
      print "\t".join(header_list)
      for sample in phenoDict_SRMK:
        print "\t".join(phenoDict_SRMK[sample])
    print >> sys.stderr, "Finished printing results\n"

  ############ self report only #####################
  elif args.proxy == "SR":  # self report only
    print >> sys.stderr, "Assigning proxy-cases based on self report (-x SR)"

    phenoDict_SR = proxy_via_selfreport(phenoDict, args.conservControl, cp, cm, cf, cs)  # self report
    print >> sys.stderr, "Finished assigning proxy-cases based on self report\n"

    if args.number is not None:
      count_relatives(args.number,kinDict,phenoDict_SR, totalCol)
      print >> sys.stderr, "Finished counting relatives of each sample who are proxy-cases or cases\n"

    if args.output == "B":
      BOLT_print(phenoDict_SRl,totalCol)
    elif args.output == "P":
      print >> sys.stderr, "This functionality not available yet\n"
    else:
      header_list=header.split("\t")
      header_list.append("F") #add new column label to header
      print "\t".join(header_list)
      for sample in phenoDict_SR:
        print "\t".join(phenoDict_SR[sample])
    print >> sys.stderr, "Finished printing results\n"

  ############# self report plus kinship ##############
  elif args.proxy == "SPK":  # self report plus kinship
    print >> sys.stderr, "Assigning proxy-cases based on self report plus kinship (-x SPK)"

    phenoDict_SR = proxy_via_selfreport(phenoDict, args.conservControl, cp, cm, cf, cs)  # self report
    print >> sys.stderr, "Finished assigning proxy-cases based on self report\n"
    phenoDict_SRPK = proxy_via_selfreport_plus_kinship(phenoDict_SR, kinDict, totalCol)  # self report plus kinships
    print >> sys.stderr, "Finished refining proxy-case assignment based on kinship\n"

    if args.number is not None:
      count_relatives(args.number, kinDict, phenoDict_SRPK, totalCol)
      print >> sys.stderr, "Finished counting relatives of each sample who are proxy-cases or cases\n"

    if args.output == "B":
      BOLT_print(phenoDict_SRPK, totalCol)
    elif args.output == "P" or args.output == "S":
      print >> sys.stderr, "This functionality not available yet\n"
    else:
      header_list=header.split("\t")
      header_list.append("F") #add new column label to header
      print "\t".join(header_list)
      for sample in phenoDict_SRPK:
        print "\t".join(phenoDict_SRPK[sample])
    print >> sys.stderr, "Finished printing results\n"

        ############### kinship only ####################
  elif args.proxy == "K":  # kinship only
    print >> sys.stderr, "Assigning proxy-cases based on kinship (-x K)"

    phenoDict_K = proxy_via_kinship(phenoDict, kinDict, args.conservControl, totalCol, cp)
    print >> sys.stderr, "Finished assigning proxy-case based on kinship\n"

    if args.number is not None:
      count_relatives(args.number,kinDict,phenoDict_K,totalCol)
      print >> sys.stderr, "Finished counting relatives of each sample who are proxy-cases or cases\n"
    
    if args.output == "B":
      BOLT_print(phenoDict_K,totalCol)
    elif args.output == "P":
      print >> sys.stderr, "This functionality not available yet\n"
    else:
      header_list=header.split("\t")
      header_list.append("F") #add new column label to header
      print "\t".join(header_list)
      for sample in phenoDict_K:
        print "\t".join(phenoDict_K[sample])
    print >> sys.stderr, "Finished printing results\n"

  ################# all proxy case designation options calculated #####################

  elif args.proxy == "A":  # all
    print >> sys.stderr, "Assigning proxy-cases based on four types of logic (-x A)"
    phenoDict_SR = proxy_via_selfreport(phenoDict, args.conservControl, cp, cm, cf, cs)  # self report
    print >> sys.stderr, "Finished assigning proxy-cases based on self report\n"
    phenoDict_SRMK = proxy_via_selfreport_minus_kinship(phenoDict_SR, kinDict, totalCol)  # self report minus kinship
    print >> sys.stderr, "Finished refining proxy-case assignment based on kinship\n"
    phenoDict_SRPK = proxy_via_selfreport_plus_kinship(phenoDict_SR, kinDict, totalCol)  # self report plus kinships
    print >> sys.stderr, "Finished refining proxy-case assignment based on kinship\n"
    phenoDict_K = proxy_via_kinship(phenoDict, kinDict, args.conservControl, totalCol, cp) #self report using only kinship
    print >> sys.stderr, "Finished assigning proxy-cases based on kinship\n"

    if args.number is not None:
      print >> sys.stderr, "This functionality is not available when -x A is invoked\n"

    if args.output == "B" or args.output == "P":
      print >> sys.stderr, "This functionality not possible when -x A is invoked\n"
    else:
      for sample in phenoDict:
        print "\t".join([sample, phenoDict_SR[sample][totalCol], phenoDict_SRMK[sample][totalCol], phenoDict_SRPK[sample][totalCol], phenoDict_K[sample][totalCol]])
    print >> sys.stderr, "Finished printing results\n"

  else:
    print >> sys.stderr, "Option for kinship designation not correct. Please use iether SMK, SR, SPK, K, or A.\n"


# call main
if __name__ == "__main__":
  main()

