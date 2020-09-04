#!/bin/bash

# Run GAT for DS pan-tissue analysis
# By Ben Laufer

echo "Sorting bed files"

sort -k1,1 -k2,2n DMRs.bed > DMRs_sorted.bed
sort -k1,1 -k2,2n background.bed > background_sorted.bed

export annotations="/share/lasallelab/Ben/DS_DBS/DMRs/GATtissue/annotations"
sort -k1,1 -k2,2n ${annotations}/FC.bed > FC_sorted.bed
sort -k1,1 -k2,2n ${annotations}/Cblm.bed > Cblm_sorted.bed 
sort -k1,1 -k2,2n ${annotations}/Glia.bed > Glia_sorted.bed 
sort -k1,1 -k2,2n ${annotations}/Neurons.bed > Neurons_sorted.bed 
sort -k1,1 -k2,2n ${annotations}/T_cells.bed > T_cells_sorted.bed 
sort -k1,1 -k2,2n ${annotations}/Fetal.bed > Fetal_sorted.bed 
sort -k1,1 -k2,2n ${annotations}/FetalCortexDMRs.bed > FetalCortexDMRs_sorted.bed
sort -k1,1 -k2,2n ${annotations}/Buccal.bed > Buccal_sorted.bed
sort -k1,1 -k2,2n ${annotations}/Blood.bed > Blood_sorted.bed
sort -k1,1 -k2,2n ${annotations}/PlacentaDMRs.bed > PlacentaDMRs_sorted.bed
sort -k1,1 -k2,2n ${annotations}/FetalCortex.bed > FetalCortex_sorted.bed
sort -k1,1 -k2,2n ${annotations}/NeonatalBlood.bed > NeonatalBlood_sorted.bed
sort -k1,1 -k2,2n ${annotations}/FC_WGBS_DMRs.bed > FC_WGBS_DMRs_sorted.bed
sort -k1,1 -k2,2n ${annotations}/neuraliPSC.bed > neuraliPSC_sorted.bed
sort -k1,1 -k2,2n ${annotations}/HeartCpGs.bed > HeartCpGs_sorted.bed
sort -k1,1 -k2,2n ${annotations}/hg38isochores.bed > hg38isochores_sorted.bed
echo "Done sorting bed files"

echo
echo "Testing for DS tissue enrichments"

module load gat/1.3.4 

call="gat-run.py \
--segments=DMRs.bed \
--with-segment-tracks \
--annotations=FC_sorted.bed \
--annotations=Cblm_sorted.bed \
--annotations=Glia_sorted.bed \
--annotations=Neurons_sorted.bed \
--annotations=T_cells_sorted.bed \
--annotations=Fetal_sorted.bed \
--annotations=FetalCortex_sorted.bed \
--annotations=Buccal_sorted.bed \
--annotations=Blood_sorted.bed \
--annotations=PlacentaDMRs_sorted.bed \
--annotations=NeonatalBlood_sorted.bed \
--annotations=FC_WGBS_DMRs_sorted.bed \
--annotations=neuraliPSC_sorted.bed \
--annotations=HeartCpGs_sorted.bed \
--workspace=background_sorted.bed \
--isochore-file=hg38isochores_sorted.bed \
--counter=segment-overlap \
--num-samples=10000 \
--num-threads=20 \
--log=DS_tissues_hyper_hypo.log \
> DS_tissues_hyper_hypo_results.tsv"

echo $call
eval $call

echo "Done testing for DS tissue enrichments"
echo

echo "Removing temporary files"
rm *_sorted.bed
echo "Done removing temporary files"

#########
# Plots #
#########

echo
echo "Plotting DS tissue enrichment results"

module load R

call="Rscript \
--vanilla \
/share/lasallelab/Ben/DS_DBS/DMRs/GATtissue/hyper_hypo/DS_GATplots.R"

echo $call
eval $call
echo
echo "Done plotting DS tissue enrichment results"

echo
echo "Script complete"
