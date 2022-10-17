#!/bin/sh
flag=0
avwcreatehd header.xml temp
../avwcreatehd++ header.xml temp2
cmp temp.nii.gz temp2.nii.gz
j=$?
echo $j
if [ $j -ne 0 ]
then
echo "Problem found in xml header create in avwcreatehd++"
flag=1
fi 

rm temp*

avwcreatehd 10 10 10 10 2 2 2 1 0 0 0 4 temp
../avwcreatehd++ 10 10 10 10 2 2 2 1 0 0 0 4 temp2
cmp temp.nii.gz temp2.nii.gz
j=$?
echo $j
if [ $j -ne 0 ]
then
echo "Problem found in manual header create in avecreatehd++"
flag=1
fi 

rm temp*

echo ""
if [ $flag -ne 0 ]
then
echo "non-zero comparison; possible problem with avwcreatehd++"
else 
echo "All comparisons zero. No problems detected with avwcreatehd++"
fi