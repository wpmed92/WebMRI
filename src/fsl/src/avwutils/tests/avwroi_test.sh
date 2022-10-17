#!/bin/sh
flag_maths=0
flag_bit=0
#test 3D volume
avwroi $FSLTESTDIR/common/avg152T1_brain temp1 15 20 15 20 15 20
../avwroi++  $FSLTESTDIR/common/avg152T1_brain temp2 15 20 15 20 15 20
cmp temp1.nii.gz temp2.nii.gz
output2=$?
if [ $output2 -ne 0 ]
then
flag_bit=1
echo "3D roi failed bit-comparison test in avwroi++"
fi
avwmaths temp1 -sub temp2 temp3
output1=`avwstats temp3 -m`
echo "Result of 3D roi test:" $output1 $output2
if [ $output1 != 0.000000 ];
then
flag_maths=1
echo "3D roi failed math test in avwroi++";
fi

rm temp*

#test 4D timeseries
avwroi $FSLTESTDIR/common/filtered_func_data.nii.gz temp1 5 10
../avwroi++ $FSLTESTDIR/common/filtered_func_data.nii.gz temp2 5 10
cmp temp1.nii.gz temp2.nii.gz
output2=$?
if [ $output2 -ne 0 ]
then
flag_bit=1
echo "4D timeseries failed bit-comparison test in avwroi++"
fi
avwmaths temp2 -Tmean temp3 
avwmaths temp1 -Tmean temp2
avwmaths temp3 -sub temp2 temp1
output1=`avwstats temp1 -m`
echo "Result of 4D timeseries test:" $output1 $output2
if [ $output1  != 0.000000 ];
then
echo "4D timeseries failed math test in avwroi++"
flag_maths=1
fi

rm temp*

#test 4D time and volume
avwroi $FSLTESTDIR/common/filtered_func_data.nii.gz temp1 15 20 15 20 5 10 5 10
../avwroi++ $FSLTESTDIR/common/filtered_func_data.nii.gz temp2 15 20 15 20 5 10 5 10
cmp temp1.nii.gz temp2.nii.gz
output2=$?
if [ $output2 -ne 0 ]
then
flag_bit=1
echo "4D timeseries and volume failed bit-comparison test in avwroi++"
fi
avwmaths temp2 -Tmean temp3 
avwmaths temp1 -Tmean temp2
avwmaths temp3 -sub temp2 temp1
output1=`avwstats temp1 -m`
echo "Result of 4D time and volume test:" $output1 $output2
if [ $output1  != 0.000000 ]
then
echo "4D timeseries and volume failed math test in avwroi++"
flag_maths=1
fi

rm temp*

echo ""
if [ $flag_maths -ne 0 -o $flag_bit -ne 0 ]
then
echo "non-zero comparison; possible problem with avwroi++"
else 
echo "All comparisons zero. No problems detected with avwroi++"
fi