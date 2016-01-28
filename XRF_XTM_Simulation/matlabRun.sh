#bin!bash
set -m 
#### When run GenerateXRFspectrum, remeber to change 'thetan=40;'
matlabFile=$1
outputFile=$2
clusterName=$(hostname)
matlab -nodesktop -nodisplay <$matlabFile &> $outputFile &
pid=$!
message=$(echo "cluster name:" $clusterName "output:" $outputFile "job id:" $pid)
mail -s 'matlab start' wendydi@mcs.anl.gov <<< $message
fg %1
mail -s 'matlab done' wendydi@mcs.anl.gov <<< $message
