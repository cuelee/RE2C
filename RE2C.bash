# !usr/bin/bash -u
# created by Cue Hyunkyu Lee
# Feb 7 2017

# ABSOLUTE_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"

# set cur directory
curdir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
parsed=""
#unparsed=""
inputf="NA"
outputf="$curdir/out"
threads="1"
corf="NA"

# load the versioninfo
cat $curdir/data/versioninfo

# read argument and parsing them
while [[ $# -gt 0 ]]
do 
key="$1"

case $key in
	--i|--input)
	inputf="$(cd "$(dirname "$2")" && pwd)/$(basename "$2")"
	parsed="$parsed input"
	shift # past argument
	;;
	--o|--output)
	outputf="$(cd "$(dirname "$2")" && pwd)/$(basename "$2")"
	parsed="$parsed output"
	shift # past argument
	;;
	--c|--cor)
	corf="$(cd "$(dirname "$2")" && pwd)/$(basename "$2")"
	parsed="$parsed cor"
	shift # past argument
	;;
	--threads)
	parsed="$parsed --threads"
	threads="$2"
	shift # past argument
	;;
	--h|--help)
	cat $curdir/data/help
	exit 0
	;;
	*)
	shift
	;;
esac # end of case

shift # past argument  
done # end of 1st while

if [ $inputf = "NA" ]
then
	echo "failed to find an input file"
	exit 0
fi

## print argument parsing information
echo " parsed arguments:$parsed"
#echo "unparsed arguments:$unparsed"
echo " ## input arguments ## "
echo " WD:	$curdir"
echo " INPUT:	$inputf"
echo " OUT:	$outputf"
echo " COR:	$corf"
echo " CORES:	$threads"


## create log file
logf="$outputf.log"
outf="$outputf.txt"

touch $outf
cat $curdir/data/versioninfo > $logf
echo " ## input arguments ## " >> $logf 
echo " WD:	$curdir" >> $logf
echo " INPUT:	$inputf" >> $logf
echo " OUTPUT:	$outputf" >> $logf
echo " COR:	$corf" >> $logf
echo " CORES:	$threads" >> $logf


## runAnalysis!!
#echo "Rscript $curdir/R/runAnalysis.R \"$curdir\" \"$inputf\" \"$outputf\" \"$corf\""
Rscript $curdir/R/runAnalysis.R $curdir $inputf $outputf $corf $logf $threads

