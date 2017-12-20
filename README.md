# proxyPower
Scripts used to analyze proxy-cases in population-based genome-wide association studies.

## Introduction

[Liu et al, 2017](https://github.com/bnwolford/proxyPower/blob/master/README.md) introduced the concept of GWAS by proxy (GWAX): performing case-control genetic association studies using unaffected first-degree relatives of cases (proxy-cases) in the (near) absence of true cases. This was motivated by large population-based biobanks, such as the UK Biobank, with a largely young and healthy population where cases for diseases such as Alzheimer's are rare, but relatives of Alzheimer's cases are present. Proxy-cases can be identified via self report of family history of disease in survey questionnaires that are frequently a part of biobank enrollment. Previous work in this area also includes the kin-cohort method ([Joshi et al, 2016](https://www.nature.com/articles/ncomms11174) and [Wacholder et al, 1998](https://academic.oup.com/aje/article/148/7/623/148336)).

The Willer group recognizes that some biobanks, such as the Nord-Trøndelag Health Study (HUNT) and the latest release of UK Biobank, will have cases and proxy-cases for diseases (e.g. type 2 diabetes or myocardial infarction). We have demonstrated the utility of these proxy-cases and developed tools for using these samples in multiple models to study the genetics of complex diseases. The scripts that will be made available here are useful for 

1. Identifying proxy-cases in your cohort
- proxyCaseAssign1dr.py 
- proxyCaseAssignAffRel.py
- proxyModel.py

2. Running statistical methods that model proxy-cases

These are beta versions and are still under development. More scripts are forthcoming.

