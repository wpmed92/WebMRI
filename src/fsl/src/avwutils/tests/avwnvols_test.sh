#!/bin/sh
#test 3D volume
flag=0
output1=`../avwnvols++ $FSLTESTDIR/common/avg152T1_brain.img` 
output2=`avwnvols $FSLTESTDIR/common/avg152T1_brain.img` 
echo "Result of 3D nvols test:" $output1 $output2
if [ $output1 != $output2 ];
then
echo "Possible problem in avwnvols++"
flag=1
fi



#test 4D timeseries
output1=`../avwnvols++ $FSLTESTDIR/common/filtered_func_data.nii.gz` 
output2=`avwnvols $FSLTESTDIR/common/filtered_func_data.nii.gz` 
echo "Result of 4D 100 frame nvols test:" $output1 $output2
if [ X$output1  != X$output2 ];
then
echo "Possible problem in avwnvols++"
flag=1
fi

echo ""
if [ $flag -ne 0 ]
then
echo "volume count discrepancy; possible problem with avwnvols++"
else 
echo "Counts are 1 and 100. No problems detected with avwnvols++"
fi