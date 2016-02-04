#bin!bash
set -m 
#### When run GenerateXRFspectrum, remeber to change 'thetan=40;'
matlabFile=$1
outputFile=$2
clusterName=$(hostname)
#  matlab -nodesktop -nodisplay <optXTM_XRF_diff.m &> opt05single10phantom_8theta_diff.out &
#  matlab -nodesktop -nodisplay <MG_simple.m &> opt23_10singlePhantom_6theta_XRF.out &
 matlab -nodesktop -nodisplay <$matlabFile &> $outputFile &
#  matlab -nodesktop -nodisplay <Conv_Test.m &> opt24_testRank.out &
# matlab -nodesktop -nodisplay <GenerateXRFspectrum.m &> opt26_xrf.out &
pid=$!
message=$(echo "cluster name:" $clusterName "output:" $outputFile "job id:" $pid)
mail -s 'matlab start' wendydi@mcs.anl.gov <<< $message
fg %1
mail -s 'matlab done' wendydi@mcs.anl.gov <<< $message
# matlab -nodesktop -nodisplay <optXTM_XRF_alternate.m &> opt05single10phantom_8theta_alter.out &
