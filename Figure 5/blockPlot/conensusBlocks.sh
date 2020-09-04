# consensusBlocks
# Create a consensus block profile from multiple comparisons
# By Ben Laufer

module load bedtools

cd /share/lasallelab/Ben/DS_DBS/DMRs

mkdir cytosine_reports/WGCNA/

echo Blocks
call="cat \
DSvsTD/Production/Blocks/blocks.bed \
DSvsDD/Production/Blocks/blocks.bed | \
sort -k1,1 -k2,2n | \
bedtools merge \
> cytosine_reports/WGCNA/ConsensusBlocks.bed"

echo $call
eval $call

echo
wc -l cytosine_reports/WGCNA/ConsensusBlocks.bed
echo
echo "Done"
echo

# Note: There were no blocks detected for DDvsTD
