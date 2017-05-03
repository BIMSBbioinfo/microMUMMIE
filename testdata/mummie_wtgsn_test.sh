#!/bin/bash
#$ -cwd
#$ -l mem_free=12G,h_vmem=20G
#$ -e error.txt
#$ -o out.txt


MUMMIE_PATH="/data/ohler/Software/microMUMMIE_new/wtgsn/" 
WD=$(pwd)	#your workig directory; has to be outside of MUMMIE_PATH
clu=$(ls *.clusters) #works only if there is just one; and it has to be PARalyzer output
PREFIX="${clu%.*}" #This determines the library prefix
TOP_miRNA=100
TwoBit="/data/ohler/harm/Genomes/drosophila/BDGP6_Ensembl_release81/TwoBit/Drosophila_melanogaster.BDGP6.dna.toplevel.2bit"
Vit_Post_toggle=0 	#0 means run viterby; 1 means run posterior


export PATH=${PATH}:${MUMMIE_PATH}
export PERL5LIB="$PERL5LIB:${MUMMIE_PATH}"
export PERL5LIB="$PERL5LIB:/home/hwessels/.guix-profile/lib/perl5/site_perl/5.24.0"
export PATH="$HOME/.guix-profile/bin:$PATH" 





##########################################################
#STEP 1
#Pre processing
#
#determines the top N number of miRNA to use from the AGO PARCLIP data
#
##########################################################


if [ ! -d ${WD}/${PREFIX}.$TOP_miRNA.1 ]
then

	#paths need to be set in path_to/preMicroMUMMIE.sh
    bash ${MUMMIE_PATH}pre-process/preMicroMUMMIE.sh ${PREFIX}.readcsv $TOP_miRNA 1
fi





##########################################################
#STEP 2
#Run microMUMMIE without targetscan
#
##########################################################

#	microMUMMIE_targetscan.pl <mature-miRNAs.txt> <mature.fa> <genome.2bit> <paralyzer-output-dir> <library-name> <out.gff> <posterior-decoding:0/1> <UTRs.txt> <Maf file> <Path of directory where you want to write all analysis in the end>\n

perl ${MUMMIE_PATH}microMUMMIE_targetscan.pl ${WD}/${PREFIX}.$TOP_miRNA.1/${PREFIX}.$TOP_miRNA.miRNA.txt ${WD}/${PREFIX}.$TOP_miRNA.1/${PREFIX}.$TOP_miRNA.miRNA.fa ${MUMMIE_PATH}Genome.2bit $WD $PREFIX ${PREFIX}.merged.gff $Vit_Post_toggle ${MUMMIE_PATH}orderedtc_scUTRs.txt ${MUMMIE_PATH}mchr_orderedtc_scUTRs.maf $WD

