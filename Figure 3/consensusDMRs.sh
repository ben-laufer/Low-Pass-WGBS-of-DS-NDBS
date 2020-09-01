# Consensus DMRs for consensusPCA.R
# By Ben Laufer

module load bedtools

cd /share/lasallelab/Ben/DS_DBS/DMRs

mkdir cytosine_reports/consensus/

echo DMRs
call="cat \
DSvsTD/Production/DMRs/DMRs.bed \
DSvsDD/Production/DMRs/DMRs.bed \
DDvsTD/Production/DMRs/DMRs.bed | \
sort -k1,1 -k2,2n | \
bedtools merge \
> cytosine_reports/consensus/ConsensusDMRs.bed"

echo $call
eval $call

echo
wc -l cytosine_reports/consensus/ConsensusDMRs.bed
echo
echo "Done"
echo
