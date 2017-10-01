humann_2_wham.py README

This script is intended to transform output gene tables from the HUMAnN2 pipeline into a wham-friendly input.
In order to execute this script successfully you must have access to python/2.7 or higher.

Execution:

python humann_2_wham.py humann2_input.tsv > wham_output.tsv

** please be sure humann_2_wham.py and the desired humann2_input is located in your working directory as you execute the script.

A sample HUMANn2 gene table will look like this...

# Gene Family	Sample1_Abundance-RPKs	Sample2_Abundance-RPKs	Sample3_Abundance-RPKs	Sample4_Abundance-RPKs
UniRef90_A0A008QTR8: NO_NAME	9.06378e-07	0	0	0
UniRef90_A0A008QTR8: NO_NAME|g__Enterococcus.s__Enterococcus_faecalis	9.06378e-07	0	0	0
UniRef90_C5EH42: Putative_stage_III_sporulation_protein_AB	1.03e-08	0	0	8.7611e-06
UniRef90_C5EH42: Putative_stage_III_sporulation_protein_AB|g__Clostridiales_noname.s__Clostridiales_bacterium_1_7_47FAA	0	0	0	8.7611e-06
UniRef90_C5EH42: Putative_stage_III_sporulation_protein_AB|unclassified	1.03e-08	0	0	0

** Notice the header line MUST start with a '#' in order for this script to identify the input as a humann2 file
** HUMANn2 organizes each gene family with a header line containing the total abundance count for all identified taxa followed by the individual taxa contributions in each subsequent line.
** This file MUST be tab delimited

The result of executing this script will look like this...

Acc	Gene_Family	Species	Sample1_Abundance-RPKs	Sample2_Abundance-RPKs	Sample3_Abundance-RPKs	Sample4_Abundance-RPKs
A0A008QTR8	NO_NAME	g__Enterococcus.s__Enterococcus_faecalis	9.06378e-07	0	0	0
C5EH42	Putative_stage_III_sporulation_protein_AB	g__Clostridiales_noname.s__Clostridiales_bacterium_1_7_47FAA	0	0	0	8.7611e-06
C5EH42	Putative_stage_III_sporulation_protein_AB	unclassified	1.03e-08	0	0	0

** Notice the header lines are removed and only the relevant lines containing abundance and species informationa are retained. This output is now ready to be uploaded to wham!

If there are any questions about this README, the humann_2_wham.py script or for general guidance in using wham! please contact the developers.
Contact information availble at our github page, https://github.com/ruggleslab/jukebox

Wham! is located at https://ruggleslab.shinyapps.io/wham_v1/