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
