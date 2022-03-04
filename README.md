# PredCRP

Predicting the regulatory role of CRP transcription factor in Escherichia coli.
This work uses an optimal feature selection method to identify 12 informative features of CRP-binding sites in cooperation with a support vector machine.
PredCRP achieved training and test accuracy of 0.98 and 0.93, respectively. This work screened and identified 23 previously unobserved regulatory interactions in Escherichia coli. PredCRP predicted the regulatory roles of CRP acting on the 23 sites and achieved test accuracy of 0.96 according to quantitative PCR validation.

Setup and Data Format
============================
On Unix systems, type `make` to build the `svm-scale` and `svm-predict`programs.  
The format of input csv file has 4 columns:  
  -  CRPBS: is a string indicating the CRP binding sites sequence from 22 to 42 base pairs by respectively adding 10 base pairs of flanking nucleotides to the regions of upstream and downstream of a CRP-binding site for considering interactions within a cis-regulatory region.  
  -  Distance of Center Position of CRPBS to TSS:  is a floating point number indicating the distance of center position of CRPBS to the transcription start site.  
  -  Transcription Unit: Transcription unit regulated by the CRP.  
  -  Regulatory Role: Gene expression effect caused by the CRP bound to the CRPBS (+ activation, - repression, ? unknown). 
  
Each line contains an instance and is ended by a '\n' character. 

Usage of PredCRP.pl
==========================
```shell
Usage: perl PrecCRP.pl [Options]  
Options:  
	-i            FILE: input CRP binding site information.  
	-svmscale     pathname: set svm-scale executable path and name (Default: svm-scale is in the same folder).  
	-svmpredict   pathname: set svm-predict executable path and name (Default: svm-scale is in the same folder).  
	-model        pathname: set PredCRP_model path and name (Default: PredCRP_model is in the same folder).  
	-h, -help
```
Usage of feature_extraction_BindingSites.pl 
===========================================
```shell
Usage: perl feature_extraction_BindingSites.pl [Options]  
Options:  
	-input		[FILE] The BindingSitesSet.txt download from RegulonDB.  
	-TF		[STR]	The interested TF (Ex: CRP).  
	-evidence	[No]	Evidence level (0:Weak, 1:Strong, 2:Both).  
	-length		[No]	10bp+BindingSites+10bp (Ex: The length of CRP is 42 (10 + 22 + 10) ).  
	-h		Show the usage.  
```

Running 23 weak-evidence data
==============================
-  step1:  Build executable programs `svm-scale` and `svm-predict`.  
```shell
$ make
```
-  step2:  The 23 weak-evidence data is in the data folder. Query input data by PredCRP.pl:  
```sh
$ perl PredCRP.pl -i data/CRPBS_23weak.csv
```
- step3: To see the prediction results  
```sh
$ cd predict_result/
$ less -S CRPBS_23weak_PredictResult.csv
```
  
Authors
=======
- Ming-Ju Tsai: milutsai.bi98g@g2.nctu.edu.tw
- Shinn-Ying Ho:syho@nctu.edu.tw

Citing PredCRP
==============
PredCRP: predicting and analysing the regulatory roles of CRP from its binding sites in Escherichia coli. Scientific Reports volume 8, Article number: 951 (2018) [PMID:[2934372](https://www.ncbi.nlm.nih.gov/pubmed/29343727)]
