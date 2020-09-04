# homer
# Parallel HOMER known transcription factor motif analyses for multiple DMRichR comparisons
# By Ben Laufer

module load homer

homerDMR(){
	i=${1}
	
	echo
	cd /share/lasallelab/Ben/DS_DBS/${i}/Production/Extra/HOMER
	echo "Running HOMER in ${PWD}"
	
    echo 
    echo "Testing hypermethylated and hypomethylated DMRs from ${i}"
	mkdir both
	cp ../../DMRs/DMRs.bed .
	
	call="findMotifsGenome.pl \
	DMRs.bed \
	hg38 \
	both/ \
	-bg background.bed \
	-cpg \
	-size given \
	-p 10 \
	-nomotif"

	echo $call
	eval $call
    
    echo 
	echo "Testing hypermethylated DMRs from ${i}"
	mkdir hyper

	call="findMotifsGenome.pl \
	DMRs_hyper.bed \
	hg38 \
	hyper/ \
	-bg background.bed \
	-cpg \
	-size given \
	-p 10 \
	-nomotif"

	echo $call
	eval $call

    echo 
    echo "Testing hypomethylated DMRs from ${i}"
	mkdir hypo

	call="findMotifsGenome.pl \
	DMRs_hypo.bed \
	hg38 \
	hypo/ \
	-bg background.bed \
	-cpg \
	-size given \
	-p 10 \
	-nomotif"

	echo $call
	eval $call
}
export -f homerDMR

parallel --jobs 3 --dry-run --will-cite  "homerDMR {}" ::: DSvsTD DSvsDD DDvsTD

parallel --jobs 3 --verbose --will-cite  "homerDMR {}" ::: DSvsTD DSvsDD DDvsTD
