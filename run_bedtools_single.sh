#!/bin/bash
patients_index=(1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62);

	for i in `seq 0 ${#patients_index[@]}`;
		do
	 	printf "%s\n" "$i";
		bedtools intersect -a ~/FinalCodeCNVs/input_bedtools_files/Benign/benign_pcnv_bed_patient_${patients_index[$i]}.bed -b ~/FinalCodeCNVs/input_bedtools_files/Benign/benign_dnacopy_bed_patient_${patients_index[$i]}.bed -wao > ~/FinalCodeCNVs/output_bedtools_files/Benign/benign_intersect_patient_wao_${patients_index[$i]}.txt;
	done

	for i in `seq 0 ${#patients_index[@]}`;
		do
	 	printf "%s\n" "$i";
		bedtools intersect -a ~/FinalCodeCNVs/input_bedtools_files/Germline/germline_pcnv_bed_patient_${patients_index[$i]}.bed -b ~/FinalCodeCNVs/input_bedtools_files/Germline/germline_dnacopy_bed_patient_${patients_index[$i]}.bed -wao > ~/FinalCodeCNVs/output_bedtools_files/Germline/germline_intersect_patient_wao_${patients_index[$i]}.txt;
	done


	for i in `seq 0 ${#patients_index[@]}`;
		do
	 	printf "%s\n" "$i";
		bedtools intersect -a ~/FinalCodeCNVs/input_bedtools_files/Tumour/tumour_pcnv_bed_patient_${patients_index[$i]}.bed -b ~/FinalCodeCNVs/input_bedtools_files/Tumour/tumour_dnacopy_bed_patient_${patients_index[$i]}.bed -wao > ~/FinalCodeCNVs/output_bedtools_files/Tumour/tumour_intersect_patient_wao_${patients_index[$i]}.txt;
	done







	for i in `seq 0 ${#patients_index[@]}`;
		do
	 	printf "%s\n" "$i";
		bedtools intersect -a ~/FinalCodeCNVs/input_bedtools_files/Benign/benign_pcnv_bed_patient_${patients_index[$i]}.bed -b ~/FinalCodeCNVs/input_bedtools_files/Germline/germline_pcnv_bed_patient_${patients_index[$i]}.bed -wao > ~/FinalCodeCNVs/output_bedtools_files/Benign_Germline/ben_germ_intersect_pcnv2_patient_wao_${patients_index[$i]}.txt;
		bedtools intersect -a ~/FinalCodeCNVs/input_bedtools_files/Benign/benign_pcnv_bed_patient_${patients_index[$i]}.bed -b ~/FinalCodeCNVs/input_bedtools_files/Germline/germline_dnacopy_bed_patient_${patients_index[$i]}.bed -wao > ~/FinalCodeCNVs/output_bedtools_files/Benign_Germline/ben_germ_intersect_pcnvDNAc_patient_wao_${patients_index[$i]}.txt;
	done




	for i in `seq 0 ${#patients_index[@]}`;
		do
	 	printf "%s\n" "$i";
		bedtools intersect -a ~/FinalCodeCNVs/input_bedtools_files/Benign/benign_pcnv_bed_patient_${patients_index[$i]}.bed -b ~/FinalCodeCNVs/input_bedtools_files/Tumour/tumour_pcnv_bed_patient_${patients_index[$i]}.bed -wao > ~/FinalCodeCNVs/output_bedtools_files/Benign_Tumour/ben_tum_intersect_pcnv2_patient_wao_${patients_index[$i]}.txt;
		bedtools intersect -a ~/FinalCodeCNVs/input_bedtools_files/Benign/benign_pcnv_bed_patient_${patients_index[$i]}.bed -b ~/FinalCodeCNVs/input_bedtools_files/Tumour/tumour_dnacopy_bed_patient_${patients_index[$i]}.bed -wao > ~/FinalCodeCNVs/output_bedtools_files/Benign_Tumour/ben_tum_intersect_pcnvDNAc_patient_wao_${patients_index[$i]}.txt;
	done





	for i in `seq 0 ${#patients_index[@]}`;
		do
	 	printf "%s\n" "$i";
		bedtools intersect -a ~/FinalCodeCNVs/input_bedtools_files/Tumour/tumour_pcnv_bed_patient_${patients_index[$i]}.bed -b ~/FinalCodeCNVs/input_bedtools_files/Germline/germline_pcnv_bed_patient_${patients_index[$i]}.bed -wao > ~/FinalCodeCNVs/output_bedtools_files/Tumour_Germline/tum_germ_intersect_pcnv2_patient_wao_${patients_index[$i]}.txt;
		bedtools intersect -a ~/FinalCodeCNVs/input_bedtools_files/Tumour/tumour_pcnv_bed_patient_${patients_index[$i]}.bed -b ~/FinalCodeCNVs/input_bedtools_files/Germline/germline_dnacopy_bed_patient_${patients_index[$i]}.bed -wao > ~/FinalCodeCNVs/output_bedtools_files/Tumour_Germline/tum_germ_intersect_pcnvDNAc_patient_wao_${patients_index[$i]}.txt;
	done







