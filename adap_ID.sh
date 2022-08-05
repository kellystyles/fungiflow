#! /bin/bash

"""
This will find adapters in you sequence.fq.gz file. Once found you can run trimmomattic or bbduk to remove them
the scrip will write the output to stdout, as well as a defined output file that can be the input to either of these progs!

Usage: ./adap_ID.sh filename.fq.gz output.fa

Written by Matt Storey, edited by Kelly Styles.
"""

#########################################################


#CLI arguements input .fq.gz file to be analysed for adapters
FILE="$1"
ADAP_FA="$2"
echo $FILE $ADAP_FA

# array of adapters (from trimmomatic etc)
declare -A ADAP=([>Reverse_adapter]="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG"
[>TruSeq_Universal_Adapter]="AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
[>pcr_dimer]="AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG"
[>PCR_Primers]="AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTCAAGCAGAAGACGGCATACGAGCTCTTCCGATCT"
[>TruSeq_Adapter_Index_1_6]="GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG"
[>TruSeq_Adapter_Index_2]="GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG"
[>TruSeq_Adapter_Index_3]="GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG"
[>TruSeq_Adapter_Index_4]="GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG"
[>TruSeq_Adapter_Index_5]="GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG"
[>TruSeq_Adapter_Index_6]="GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG"
[>TruSeq_Adapter_Index_7]="GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG"
[>TruSeq_Adapter_Index_8]="GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG"
[>TruSeq_Adapter_Index_9]="GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG"
[>TruSeq_Adapter_Index_10]="GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG"
[>TruSeq_Adapter_Index_11]="GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG"
[>TruSeq_Adapter_Index_12]="GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG"
[>TruSeq_Adapter_Index_13]="GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACAATCTCGTATGCCGTCTTCTGCTTG"
[>TruSeq_Adapter_Index_14]="GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTATCTCGTATGCCGTCTTCTGCTTG"
[>TruSeq_Adapter_Index_15]="GATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAGAATCTCGTATGCCGTCTTCTGCTTG"
[>TruSeq_Adapter_Index_16]="GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCCGATCTCGTATGCCGTCTTCTGCTTG"
[>TruSeq_Adapter_Index_18_7]="GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCACATCTCGTATGCCGTCTTCTGCTTG"
[>TruSeq_Adapter_Index_19]="GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACGATCTCGTATGCCGTCTTCTGCTTG"
[>TruSeq_Adapter_Index_20]="GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGGCCTTATCTCGTATGCCGTCTTCTGCTTG"
[>TruSeq_Adapter_Index_21]="GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTTTCGGAATCTCGTATGCCGTCTTCTGCTTG"
[>TruSeq_Adapter_Index_22]="GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTAATCTCGTATGCCGTCTTCTGCTTG"
[>TruSeq_Adapter_Index_23]="GATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTGGATATCTCGTATGCCGTCTTCTGCTTG"
[>TruSeq_Adapter_Index_25]="GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATATCTCGTATGCCGTCTTCTGCTTG"
[>TruSeq_Adapter_Index_27]="GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTTTATCTCGTATGCCGTCTTCTGCTTG"
[>I5_Nextera_Transposase_1]="CTGTCTCTTATACACATCTGACGCTGCCGACGA"
[>I7_Nextera_Transposase_1]="CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
[>I5_Nextera_Transposase_2]="CTGTCTCTTATACACATCTCTGATGGCGCGAGGGAGGC"
[>I7_Nextera_Transposase_2]="CTGTCTCTTATACACATCTCTGAGCGGGCTGGCAAGGC"
[>I5_Primer_Nextera_XT_and_Nextera_Enrichment_NSE_501]="GACGCTGCCGACGAGCGATCTAGTGTAGATCTCGGTGGTCGCCGTATCATT"
[>I5_Primer_Nextera_XT_and_Nextera_Enrichment_NSE_502]="GACGCTGCCGACGAATAGAGAGGTGTAGATCTCGGTGGTCGCCGTATCATT"
[>I5_Primer_Nextera_XT_and_Nextera_Enrichment_NSE_503]="GACGCTGCCGACGAAGAGGATAGTGTAGATCTCGGTGGTCGCCGTATCATT"
[>I5_Primer_Nextera_XT_and_Nextera_Enrichment_NSE_504]="GACGCTGCCGACGATCTACTCTGTGTAGATCTCGGTGGTCGCCGTATCATT"
[>I5_Primer_Nextera_XT_and_Nextera_Enrichment_NSE_505]="GACGCTGCCGACGACTCCTTACGTGTAGATCTCGGTGGTCGCCGTATCATT"
[>I5_Primer_Nextera_XT_and_Nextera_Enrichment_NSE_506]="GACGCTGCCGACGATATGCAGTGTGTAGATCTCGGTGGTCGCCGTATCATT"
[>I5_Primer_Nextera_XT_and_Nextera_Enrichment_NSE_507]="GACGCTGCCGACGATACTCCTTGTGTAGATCTCGGTGGTCGCCGTATCATT"
[>I5_Primer_Nextera_XT_and_Nextera_Enrichment_NSE_508]="GACGCTGCCGACGAAGGCTTAGGTGTAGATCTCGGTGGTCGCCGTATCATT"
[>I5_Primer_Nextera_XT_and_Nextera_Enrichment_NSE_517]="GACGCTGCCGACGATCTTACGCGTGTAGATCTCGGTGGTCGCCGTATCATT"
[>I7_Primer_Nextera_XT_and_Nextera_Enrichment_N701]="CCGAGCCCACGAGACTAAGGCGAATCTCGTATGCCGTCTTCTGCTTG"
[>I7_Primer_Nextera_XT_and_Nextera_Enrichment_N702]="CCGAGCCCACGAGACCGTACTAGATCTCGTATGCCGTCTTCTGCTTG"
[>I7_Primer_Nextera_XT_and_Nextera_Enrichment_N703]="CCGAGCCCACGAGACAGGCAGAAATCTCGTATGCCGTCTTCTGCTTG"
[>I7_Primer_Nextera_XT_and_Nextera_Enrichment_N704]="CCGAGCCCACGAGACTCCTGAGCATCTCGTATGCCGTCTTCTGCTTG"
[>I7_Primer_Nextera_XT_and_Nextera_Enrichment_N705]="CCGAGCCCACGAGACGGACTCCTATCTCGTATGCCGTCTTCTGCTTG"
[>I7_Primer_Nextera_XT_and_Nextera_Enrichment_N706]="CCGAGCCCACGAGACTAGGCATGATCTCGTATGCCGTCTTCTGCTTG"
[>I7_Primer_Nextera_XT_and_Nextera_Enrichment_N707]="CCGAGCCCACGAGACCTCTCTACATCTCGTATGCCGTCTTCTGCTTG"
[>I7_Primer_Nextera_XT_and_Nextera_Enrichment_N708]="CCGAGCCCACGAGACCAGAGAGGATCTCGTATGCCGTCTTCTGCTTG"
[>I7_Primer_Nextera_XT_and_Nextera_Enrichment_N709]="CCGAGCCCACGAGACGCTACGCTATCTCGTATGCCGTCTTCTGCTTG"
[>I7_Primer_Nextera_XT_and_Nextera_Enrichment_N710]="CCGAGCCCACGAGACCGAGGCTGATCTCGTATGCCGTCTTCTGCTTG"
[>I7_Primer_Nextera_XT_and_Nextera_Enrichment_N711]="CCGAGCCCACGAGACAAGAGGCAATCTCGTATGCCGTCTTCTGCTTG"
[>I7_Primer_Nextera_XT_and_Nextera_Enrichment_N712]="CCGAGCCCACGAGACGTAGAGGAATCTCGTATGCCGTCTTCTGCTTG"
[>I5_Primer_Nextera_XT_Index_Kit_v2_S502]="GACGCTGCCGACGAATAGAGAGGTGTAGATCTCGGTGGTCGCCGTATCATT"
[>I5_Primer_Nextera_XT_Index_Kit_v2_S503]="GACGCTGCCGACGAAGAGGATAGTGTAGATCTCGGTGGTCGCCGTATCATT"
[>I5_Primer_Nextera_XT_Index_Kit_v2_S505]="GACGCTGCCGACGACTCCTTACGTGTAGATCTCGGTGGTCGCCGTATCATT"
[>I5_Primer_Nextera_XT_Index_Kit_v2_S506]="GACGCTGCCGACGATATGCAGTGTGTAGATCTCGGTGGTCGCCGTATCATT"
[>I5_Primer_Nextera_XT_Index_Kit_v2_S507]="GACGCTGCCGACGATACTCCTTGTGTAGATCTCGGTGGTCGCCGTATCATT"
[>I5_Primer_Nextera_XT_Index_Kit_v2_S508]="GACGCTGCCGACGAAGGCTTAGGTGTAGATCTCGGTGGTCGCCGTATCATT"
[>I5_Primer_Nextera_XT_Index_Kit_v2_S510]="GACGCTGCCGACGAATTAGACGGTGTAGATCTCGGTGGTCGCCGTATCATT"
[>I5_Primer_Nextera_XT_Index_Kit_v2_S511]="GACGCTGCCGACGACGGAGAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
[>I5_Primer_Nextera_XT_Index_Kit_v2_S513]="GACGCTGCCGACGACTAGTCGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
[>I5_Primer_Nextera_XT_Index_Kit_v2_S515]="GACGCTGCCGACGAAGCTAGAAGTGTAGATCTCGGTGGTCGCCGTATCATT"
[>I5_Primer_Nextera_XT_Index_Kit_v2_S516]="GACGCTGCCGACGAACTCTAGGGTGTAGATCTCGGTGGTCGCCGTATCATT"
[>I5_Primer_Nextera_XT_Index_Kit_v2_S517]="GACGCTGCCGACGATCTTACGCGTGTAGATCTCGGTGGTCGCCGTATCATT"
[>I5_Primer_Nextera_XT_Index_Kit_v2_S518]="GACGCTGCCGACGACTTAATAGGTGTAGATCTCGGTGGTCGCCGTATCATT"
[>I5_Primer_Nextera_XT_Index_Kit_v2_S520]="GACGCTGCCGACGAATAGCCTTGTGTAGATCTCGGTGGTCGCCGTATCATT"
[>I5_Primer_Nextera_XT_Index_Kit_v2_S521]="GACGCTGCCGACGATAAGGCTCGTGTAGATCTCGGTGGTCGCCGTATCATT"
[>I5_Primer_Nextera_XT_Index_Kit_v2_S522]="GACGCTGCCGACGATCGCATAAGTGTAGATCTCGGTGGTCGCCGTATCATT"
[>I7_Primer_Nextera_XT_Index_Kit_v2_N701]="CCGAGCCCACGAGACTAAGGCGAATCTCGTATGCCGTCTTCTGCTTG"
[>I7_Primer_Nextera_XT_Index_Kit_v2_N702]="CCGAGCCCACGAGACCGTACTAGATCTCGTATGCCGTCTTCTGCTTG"
[>I7_Primer_Nextera_XT_Index_Kit_v2_N703]="CCGAGCCCACGAGACAGGCAGAAATCTCGTATGCCGTCTTCTGCTTG"
[>I7_Primer_Nextera_XT_Index_Kit_v2_N704]="CCGAGCCCACGAGACTCCTGAGCATCTCGTATGCCGTCTTCTGCTTG"
[>I7_Primer_Nextera_XT_Index_Kit_v2_N705]="CCGAGCCCACGAGACGGACTCCTATCTCGTATGCCGTCTTCTGCTTG"
[>I7_Primer_Nextera_XT_Index_Kit_v2_N706]="CCGAGCCCACGAGACTAGGCATGATCTCGTATGCCGTCTTCTGCTTG"
[>I7_Primer_Nextera_XT_Index_Kit_v2_N707]="CCGAGCCCACGAGACCTCTCTACATCTCGTATGCCGTCTTCTGCTTG"
[>I7_Primer_Nextera_XT_Index_Kit_v2_N710]="CCGAGCCCACGAGACCGAGGCTGATCTCGTATGCCGTCTTCTGCTTG"
[>I7_Primer_Nextera_XT_Index_Kit_v2_N711]="CCGAGCCCACGAGACAAGAGGCAATCTCGTATGCCGTCTTCTGCTTG"
[>I7_Primer_Nextera_XT_Index_Kit_v2_N712]="CCGAGCCCACGAGACGTAGAGGAATCTCGTATGCCGTCTTCTGCTTG"
[>I7_Primer_Nextera_XT_Index_Kit_v2_N714]="CCGAGCCCACGAGACGCTCATGAATCTCGTATGCCGTCTTCTGCTTG"
[>I7_Primer_Nextera_XT_Index_Kit_v2_N715]="CCGAGCCCACGAGACATCTCAGGATCTCGTATGCCGTCTTCTGCTTG"
[>I7_Primer_Nextera_XT_Index_Kit_v2_N716]="CCGAGCCCACGAGACACTCGCTAATCTCGTATGCCGTCTTCTGCTTG"
[>I7_Primer_Nextera_XT_Index_Kit_v2_N718]="CCGAGCCCACGAGACGGAGCTACATCTCGTATGCCGTCTTCTGCTTG"
[>I7_Primer_Nextera_XT_Index_Kit_v2_N719]="CCGAGCCCACGAGACGCGTAGTAATCTCGTATGCCGTCTTCTGCTTG"
[>I7_Primer_Nextera_XT_Index_Kit_v2_N720]="CCGAGCCCACGAGACCGGAGCCTATCTCGTATGCCGTCTTCTGCTTG"
[>I7_Primer_Nextera_XT_Index_Kit_v2_N721]="CCGAGCCCACGAGACTACGCTGCATCTCGTATGCCGTCTTCTGCTTG"
[>I7_Primer_Nextera_XT_Index_Kit_v2_N722]="CCGAGCCCACGAGACATGCGCAGATCTCGTATGCCGTCTTCTGCTTG"
[>I7_Primer_Nextera_XT_Index_Kit_v2_N723]="CCGAGCCCACGAGACTAGCGCTCATCTCGTATGCCGTCTTCTGCTTG"
[>I7_Primer_Nextera_XT_Index_Kit_v2_N724]="CCGAGCCCACGAGACACTGAGCGATCTCGTATGCCGTCTTCTGCTTG"
[>I7_Primer_Nextera_XT_Index_Kit_v2_N726]="CCGAGCCCACGAGACCCTAAGACATCTCGTATGCCGTCTTCTGCTTG"
[>I7_Primer_Nextera_XT_Index_Kit_v2_N727]="CCGAGCCCACGAGACCGATCAGTATCTCGTATGCCGTCTTCTGCTTG"
[>I7_Primer_Nextera_XT_Index_Kit_v2_N728]="CCGAGCCCACGAGACTGCAGCTAATCTCGTATGCCGTCTTCTGCTTG"
[>I7_Primer_Nextera_XT_Index_Kit_v2_N729]="CCGAGCCCACGAGACTCGACGTCATCTCGTATGCCGTCTTCTGCTTG"
[>I5_Adapter_Nextera]="CTGATGGCGCGAGGGAGGCGTGTAGATCTCGGTGGTCGCCGTATCATT"
[>I7_Adapter_Nextera_No_Barcode]="CTGAGCGGGCTGGCAAGGCAGACCGATCTCGTATGCCGTCTTCTGCTTG"
[>Nextera_LMP_Read1_External_Adapter]="GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
[>Nextera_LMP_Read2_External_Adapter]="GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
[>RNA_Adapter_RA5_part_#_15013205]="GATCGTCGGACTGTAGAACTCTGAAC"
[>RNA_Adapter_RA3_part_#_15013207]="CCTTGGCACCCGAGAATTCCA"
[>Stop_Oligo_STP_8]="CCACGGGAACGTGGTGGAATTC"
[>RNA_RT_Primer_RTP_part_#_15013981]="TGGAATTCTCGGGTGCCAAGGC"
[>RNA_PCR_Primer_RP1_part_#_15013198]="TCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT"
[>RNA_PCR_Primer_Index_1_RPI1_2,9]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_2_RPI2]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_3_RPI3]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_4_RPI4]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_5_RPI5]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_6_RPI6]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_7_RPI7]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_8_RPI8]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_9_RPI9]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_10_RPI10]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_11_RPI11]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_12_RPI12]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_13_RPI13]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACAGTCAAATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_14_RPI14]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACAGTTCCATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_15_RPI15]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATGTCAATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_16_RPI16]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCCGTCCATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_17_RPI17]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTAGAGATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_18_RPI18]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_19_RPI19]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTGAAAATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_20_RPI20]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTGGCCATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_22_RPI22]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCGTACGATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_23_RPI23]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGAGTGGATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_24_RPI24]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGGTAGCATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_25_RPI25]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACACTGATATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_26_RPI26]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATGAGCATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_27_RPI27]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATTCCTATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_28_RPI28]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCAAAAGATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_29_RPI29]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCAACTAATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_30_RPI30]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCACCGGATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_31_RPI31]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCACGATATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_32_RPI32]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCACTCAATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_33_RPI33]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCAGGCGATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_34_RPI34]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCATGGCATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_35_RPI35]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCATTTTATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_36_RPI36]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCCAACAATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_37_RPI37]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCGGAATATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_38_RPI38]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCTAGCTATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_39_RPI39]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCTATACATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_40_RPI40]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCTCAGAATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_41_RPI41]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGACGACATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_42_RPI42]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTAATCGATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_43_RPI43]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTACAGCATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_44_RPI44]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTATAATATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_45_RPI45]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTCATTCATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_46_RPI46]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTCCCGAATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_47_RPI47]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTCGAAGATCTCGTATGCCGTCTTCTGCTTG"
[>RNA_PCR_Primer_Index_48_RPI48]="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTCGGCAATCTCGTATGCCGTCTTCTGCTTG"
[>PhiX_read1_adapter]="AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTGAAA"
[>PhiX_read2_adapter]="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAA"
[>Bisulfite_R1]="AGATCGGAAGAGCACACGTCTGAAC"
[>Bisulfite_R2]="AGATCGGAAGAGCGTCGTGTAGGGA")

#loop opver the array of adapters and grep them against the input file. zcat takes .gz as imput, could set up check for file type and make allowences for uncompressed files?
for i in ${!ADAP[@]}; do
    VAL="$( zcat $FILE | head -400000 | grep "${ADAP[$i]}" | wc -l )"
    if [ $VAL -gt 100 ]
    then
        echo "$VAL"
        echo "This adapter was found:"
        echo "$i"
        echo "${ADAP[$i]}"
        echo "$i" >> $ADAP_FA
        echo "${ADAP[$i]}" >> $ADAP_FA
   fi
done
