#!/bin/sh
#test 3D volume
echo "Running awcc non-compliant voxel dimension test"
avwcc $FSLTESTDIR/common/avg152T1_brain.img $FSLTESTDIR/common/filtered_func_data.nii.gz > test
../avwcc++ $FSLTESTDIR/common/avg152T1_brain.img $FSLTESTDIR/common/filtered_func_data.nii.gz > test2
flag=0 
cmp test test2
j=$?
echo $j
if [ $j -ne 0 ]
then
echo "Problem found in this comparision"
flag=1
fi 

rm test*

avwcc $FSLTESTDIR/common/filtered_func_data.nii.gz $FSLTESTDIR/common/filtered_func_data.nii.gz > test
../avwcc++ $FSLTESTDIR/common/filtered_func_data.nii.gz $FSLTESTDIR/common/filtered_func_data.nii.gz > test2
flag=0 
cmp test test2
j=$?
echo $j
if [ $j -ne 0 ]
then
echo "Problem found in this comparision"
flag=1
fi 

rm test*

echo ""
if [ $flag -ne 0 ]
then
echo "non-zero comparison; possible problem with avwcc"
else 
echo "All comparisons zero. No problems detected with avwcc"
fi