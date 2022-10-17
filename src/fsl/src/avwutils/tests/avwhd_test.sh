#!/bin/sh
avwhd $FSLTESTDIR/common/avg152T1_brain.hdr > test
../avwhd++ $FSLTESTDIR/common/avg152T1_brain.hdr > test2
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

avwhd $FSLTESTDIR/common/filtered_func_data.nii.gz   > test
../avwhd++  $FSLTESTDIR/common/filtered_func_data.nii.gz > test2
cmp test test2
j=$?
echo $j
if [ $j -ne 0 ]
then
echo "Problem found in this comparision"
flag=1
fi 

rm test*

avwhd -x $FSLTESTDIR/common/avg152T1_brain.hdr > test
../avwhd++ -x $FSLTESTDIR/common/avg152T1_brain.hdr > test2
cmp test test2
j=$?
echo $j
if [ $j -ne 0 ]
then
echo "Problem found in this comparision"
flag=1
fi 

rm test*

avwhd -x $FSLTESTDIR/common/filtered_func_data.nii.gz   > test
../avwhd++ -x $FSLTESTDIR/common/filtered_func_data.nii.gz > test2
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
echo "non-zero comparison; possible problem with avwhd++"
else 
echo "All comparisons zero. No problems detected with avwhd++"
fi