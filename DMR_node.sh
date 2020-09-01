# Parallel DM.R runs for a high resource node
# Run from one level above the cytosine_reports folder with sample_info.xlsx files for each contrast name (i.e. sample_info_DSvsTD.xlsx)
# By Ben Laufer

# Main function

DMR(){

	contrast=$1

	mkdir ${contrast}
	cp cytosine_reports/*.gz ${contrast}
	mv sample_info_${contrast}.xlsx ${contrast}/sample_info.xlsx

	cd ${contrast}

	call="Rscript \
	--vanilla \
	/share/lasallelab/programs/DMRichR/DM.R \
	--genome hg38 \
	--coverage 1 \
	--perGroup '0.75' \
	--minCpGs 5 \
	--maxPerms 10 \
	--cutoff '0.05' \
	--testCovariate Diagnosis \
	--adjustCovariate Sex \
	--cores 20"

	echo $call
	eval $call

	rm *.gz

	echo ${contrast} contrast complete

}
export -f DMR

# Run the function using GNU parallel

module load R/3.6.3 

cd /share/lasallelab/Ben/DS_DBS/

mkdir dmrLogs

parallel --dry-run --will-cite --results dmrLogs -j 3 "DMR {}" ::: DSvsTD DSvsDD DDvsTD

parallel --verbose --will-cite --results dmrLogs -j 3 "DMR {}" ::: DSvsTD DSvsDD DDvsTD

