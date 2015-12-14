#---------------------------------- Bedtools - Installed only on the mac --------------------------------------------------#

# run_bedtools.sh - in the input_bedtools_files folder - !!! Make sure you change the patients arrays in the sh file!!!

uni_pat_bg_tg <- uni_pat_bg[uni_pat_bg%in%uni_pat_tg]

uni_pat_bg_tb <- uni_pat_bg[uni_pat_bg%in%uni_pat_tb]
# If patient x is in bg but not in tb it will still return no intersection but if I select it before it will take less time.

#--------------------------------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------------------------------#
# Read the output of bedtools
#--------------------------------------------------------------------------------------------------------------------------#

source("functions_list.R")

#BGpcnv vs BGdnacopy
create_file_with_values(uni_pat_bg,
                        "~/FinalCodeCNVs/output_bedtools_files/BG/bg_intersect_patient_wao_",
                        "~/FinalCodeCNVs/input_bedtools_files/BG/bg_dnacopy_bed_values_patient_",
                        "~/FinalCodeCNVs/input_bedtools_files/BG/bg_pcnv_bed_values_patient_",
                        "~/FinalCodeCNVs/output_bedtools_files/BGvalues/bg_intersect_with_values_patient_")
#TGpcnv vs TGdnacopy
create_file_with_values(uni_pat_tg,
                        "~/FinalCodeCNVs/output_bedtools_files/TG/tg_intersect_patient_wao_",
                        "~/FinalCodeCNVs/input_bedtools_files/TG/tg_dnacopy_bed_values_patient_",
                        "~/FinalCodeCNVs/input_bedtools_files/TG/tg_pcnv_bed_values_patient_",
                        "~/FinalCodeCNVs/output_bedtools_files/TGvalues/tg_intersect_with_values_patient_")
#TBpcnv vs TBdnacopy
create_file_with_values(uni_pat_tb,
                        "~/FinalCodeCNVs/output_bedtools_files/TB/tb_intersect_patient_wao_",
                        "~/FinalCodeCNVs/input_bedtools_files/TB/tb_dnacopy_bed_values_patient_",
                        "~/FinalCodeCNVs/input_bedtools_files/TB/tb_pcnv_bed_values_patient_",
                        "~/FinalCodeCNVs/output_bedtools_files/TBvalues/tb_intersect_with_values_patient_")
#BGpcnv vs TGpcnv
create_file_with_values(uni_pat_bg_tg,
                        "~/FinalCodeCNVs/output_bedtools_files/BG_TG/bg_tg_intersect_pcnv2_patient_wao_",
                        "~/FinalCodeCNVs/input_bedtools_files/TG/tg_pcnv_bed_values_patient_",
                        "~/FinalCodeCNVs/input_bedtools_files/BG/bg_pcnv_bed_values_patient_",
                        "~/FinalCodeCNVs/output_bedtools_files/BG_TGvalues/output_bg_tg_intersect_pcnv2_with_values_patient_")
#BGpcnv vs TBpcnv
create_file_with_values(uni_pat_bg_tb,
                        "~/FinalCodeCNVs/output_bedtools_files/BG_TB/bg_tb_intersect_pcnv2_patient_wao_",
                        "~/FinalCodeCNVs/input_bedtools_files/TB/tb_pcnv_bed_values_patient_",
                        "~/FinalCodeCNVs/input_bedtools_files/BG/bg_pcnv_bed_values_patient_",
                        "~/FinalCodeCNVs/output_bedtools_files/BG_TBvalues/output_bg_tb_intersect_pcnv2_with_values_patient_")

#BGpcnv vs TGdnacopy
create_file_with_values(uni_pat_bg_tg,
                        "~/FinalCodeCNVs/output_bedtools_files/BG_TG/bg_tg_intersect_pcnvDNAc_patient_wao_",
                        "~/FinalCodeCNVs/input_bedtools_files/TG/tg_dnacopy_bed_values_patient_",
                        "~/FinalCodeCNVs/input_bedtools_files/BG/bg_pcnv_bed_values_patient_",
                        "~/FinalCodeCNVs/output_bedtools_files/BG_TGvalues/output_bg_tg_intersect_pcnvDNAc_with_values_patient_")
#BGpcnv vs TBdnacopy
create_file_with_values(uni_pat_bg_tb,
                        "~/FinalCodeCNVs/output_bedtools_files/BG_TB/bg_tb_intersect_pcnvDNAc_patient_wao_",
                        "~/FinalCodeCNVs/input_bedtools_files/TB/tb_dnacopy_bed_values_patient_",
                        "~/FinalCodeCNVs/input_bedtools_files/BG/bg_pcnv_bed_values_patient_",
                        "~/FinalCodeCNVs/output_bedtools_files/BG_TBvalues/output_bg_tb_intersect_pcnvDNAc_with_values_patient_")


#---------------------------------------------------------------------------------------------------------------------------#
#  Analyse the txt file and build the contingency table
#---------------------------------------------------------------------------------------------------------------------------#

which_patients_pcnv(output_penncnv_bg,uni_pat_bg,"~/FinalCodeCNVs/output_bedtools_files/BG_pcnv_patients_chrs.csv")
which_patients_pcnv(output_penncnv_tg,uni_pat_tg,"~/FinalCodeCNVs/output_bedtools_files/TG_pcnv_patients_chrs.csv")
which_patients_pcnv(output_penncnv_tb,uni_pat_tb,"~/FinalCodeCNVs/output_bedtools_files/TB_pcnv_patients_chrs.csv")


setwd("~/FinalCodeCNVs")

source("functions_list.R")

create_contingency(uni_pat_bg,
                   "~/FinalCodeCNVs/output_bedtools_files/BGvalues/bg_intersect_with_values_patient_",
                   "~/FinalCodeCNVs/output_bedtools_files/BG_contingency.csv")
create_contingency(uni_pat_tg,
                   "~/FinalCodeCNVs/output_bedtools_files/TGvalues/tg_intersect_with_values_patient_",
                   "~/FinalCodeCNVs/output_bedtools_files/TG_contingency.csv")
create_contingency(uni_pat_tb,
                   "~/FinalCodeCNVs/output_bedtools_files/TBvalues/tb_intersect_with_values_patient_",
                   "~/FinalCodeCNVs/output_bedtools_files/TB_contingency.csv")

#----------------------------------------------------------------------------------------------

# Read the compare files and add the corresponding values

create_contingency_compare(uni_pat_bg_tg, 
                           "~/FinalCodeCNVs/output_bedtools_files/BG_TGvalues/output_bg_tg_intersect_pcnv2_with_values_patient_",
                           "~/FinalCodeCNVs/output_bedtools_files/BG_TGvalues/output_bg_tg_intersect_pcnvDNAc_with_values_patient_",
                           "~/FinalCodeCNVs/output_bedtools_files/contingency_bg_tg.csv", 
                           c("TG_CN=1", "TG_CN=3", "TG_CN=4","TG-Gain=1", "TG-Loss=-1", "TG-Gain=0 , TG-Loss=0"))

create_contingency_compare(uni_pat_bg_tb, 
                           "~/FinalCodeCNVs/output_bedtools_files/BG_TBvalues/output_bg_tb_intersect_pcnv2_with_values_patient_",
                           "~/FinalCodeCNVs/output_bedtools_files/BG_TBvalues/output_bg_tb_intersect_pcnvDNAc_with_values_patient_",
                           "~/FinalCodeCNVs/output_bedtools_files/contingency_bg_tb.csv", 
                           c("TB_CN=1", "TB_CN=3", "TB_CN=4","TB-Gain=1", "TB-Loss=-1", "TB-Gain=0 , TB-Loss=0"))
