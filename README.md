# trio_dataset
Copy number variation analysis of SNPs microarray from the CamCap and ICGC project.

Data

CamCap has three samples “germline, benign, tumour
ICGC has only germline (Blood) and Tumour. 
Where CC doesn’t have the blood sample, this is taken from the ICGC.

The table includes a total number of 74 patients of which 28 have all three samples for the CamCap project.

The CamCap data has already been processed using GenomeStudio. So we already have the .txt files from Neal group.
The ICGC data, instead, only have the .IDAT file so I processes them using GenomeStudio (For details see genomestudio_guide.doc inside the Output_ProstateTrio_GenomeStudio folder).
For a quicker analysis I run Genome Studio with multiple datasets (it makes no difference in terms of the results obtained but it’s much quicker than run the software to each sample separately). 
To separate each file then run the script:
FinalCodeCNVs/GenomeStudio_into_singlefiles.R 

Now I am ready to analyse the data.

#
# final_code1.R
#

Pre-processing


The file Arraystolookat.csv is a summary table, which includes the correspondences between patient samples and files, and the project.
Each patient should have at least 3 samples: germline, benign and tumour, for CamCap or germline-tumour for ICGC and benign-tumour for CC. Or I can also have all three for CC and both for ICGC. In that case you can compare the tumour signal in CC with the one in CC and see if it is the same – check reliability.

1- Discard NAs

The samples all contains in both Logr and BAF some NAs. For a better analysis I discard all NAs and put the mean of the values around (window of 3) .

2- Add 50 simulated data points

Based on the paper: “Breaking the waves: improved detection of copy number variation form microarray-based comparative genomic hybridization”, I added 50 data points before at the beginning and at the end of each chromosome.

3- LOESS

Because the data are extremely noisy, I try to correct it using loess (cit). Tis technique simply interpolates the signal and return a smoother one.

4- Eliminate the 50 simulated datapoints

5-  Runmed – window=5

At this point I have two options
-	I calculate the ratio between Benign and Germline (LogBenign – LogGermline) and analyse thi new signal
-	I leave the two signals separate and analyse them individually


