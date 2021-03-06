# proxyPower (v1.0.0)
This repository contains scripts used to analyze proxy-cases in population-based genome-wide association studies. 
:warning:These scripts are beta versions and this respository is under development.  

## Introduction

[Liu et al, 2017](https://github.com/bnwolford/proxyPower/blob/master/README.md) introduced the concept of GWAS by proxy (GWAX): performing case-control genetic association studies using unaffected first-degree relatives of cases (proxy-cases) in the (near) absence of true cases. This was motivated by large population-based biobanks, such as the UK Biobank, with a largely young and healthy population where cases for diseases such as Alzheimer's are rare, but relatives of Alzheimer's cases are present. Proxy-cases can be identified via self report of family history of disease in survey questionnaires that are frequently a part of biobank enrollment. Previous work in this area also includes the kin-cohort method ([Joshi et al, 2016](https://www.nature.com/articles/ncomms11174) and [Wacholder et al, 1998](https://academic.oup.com/aje/article/148/7/623/148336)). 

Population-based cohorts such as the Nord-Trøndelag Health Study (HUNT) and the UK Biobank have cases and proxy-cases for many diseases of interest (e.g. type 2 diabetes or myocardial infarction). We have demonstrated the utility of these proxy-cases and developed tools for using these samples in multiple models to study the genetics of complex diseases. Please view the [provided PDF](https://github.com/bnwolford/proxyPower/blob/master/proxyPower.pdf) for more information. 

## Contents

Scripts made available here are useful for:

1. Identifying proxy-cases in your cohort
- proxyCaseAssign1dr.py 
- proxyCaseAssignAffRel.py
- proxyModel.py

2. Running statistical methods that model proxy-cases

## Getting Started

In order to download proxyPower you should clone this repository.

    git clone https:://github.com/bnwolford/proxyPower.git
    cd proxyPower
 
## Support
 
 - [Tutorial](https://github.com/bnwolford/proxyPower/wiki/Tutorial)
 - Issues with proxyPower? Email bwolford@umich.edu. 
 
## License
 
This project is licensed under the terms of the MIT license.

