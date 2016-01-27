# RiboSeq

Manage ribo-seq experiment from raw data


## GSE.py

usage: GSE.py -g GSEnumber
option: -w workdir
  
  
### What is GSE number:
GSE are GEO database identification. It is associated to a list of sample
used during the study. This number can be found on the GEO NCBI database  
and is often referred on the publication.
  
### What GSE.py do:
From a GSE number, the script download the xml file containing the description of  
the experiment. After parsing the xml file GSEnumber.xml, we extract the name of  
each sample, the species, the GSM and SRX number. These latter numbers are required
to download SRA files associated.
  
### What is the output:
At the end, the user obtains a GSEnumber.xml and a GSEnumber.sample.  
+ GSEnumber.xml is the raw file of the experiment. It contains all the details, protocoles etc...  
+ GSEnumber.sample resume the samples by giving GSM, sample's names, SRR and Genome assmebly. This file is required to download SRA data.  
Output files is created in the current directory unless user specified -w option.
  
  
