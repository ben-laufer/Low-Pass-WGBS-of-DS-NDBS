# Parallel runs of Down syndrome pan-tissue analysis with GAT
# If you modify this you will also have to change the scripts and files in the other folders
# By Ben Laufer

module load R
module load gat

DSGAT(){
	i=${1}
	
	echo
	cd /share/lasallelab/Ben/DS_DBS/DMRs/${i}/Production/Extra/GAT
	echo "Running DS GAT Tissue for ${i} in ${PWD}"
	
	echo "Analyzing hypermethylated and hypomethylated separately"
	/share/lasallelab/Ben/DS_DBS/DMRs/GATtissue/hyper_hypo/DS_GAT.sh
	
	echo "Analyzing hypermethylated and hypomethylated together"
	/share/lasallelab/Ben/DS_DBS/DMRs/GATtissue/both/DS_GAT.sh
	
	echo "DS GAT Tissue has finished for ${i}"
	echo
	}
export -f DSGAT
	
parallel --dry-run --will-cite --jobs 3 "DSGAT {}" ::: DSvsTD DSvsDD DDvsTD

parallel --verbose --will-cite --jobs 3 "DSGAT {}" ::: DSvsTD DSvsDD DDvsTD
