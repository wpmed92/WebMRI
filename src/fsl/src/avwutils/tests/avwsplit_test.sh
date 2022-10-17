#!/bin/sh
#test 3D volume
avwsplit $FSLTESTDIR/common/filtered_func_data.nii.gz 
../avwsplit++ $FSLTESTDIR/common/filtered_func_data.nii.gz test
flag=0
i=0
while [ $i -lt 100 ] 
do
if [ $i -lt 10 ] 
then 
avwmaths  vol000${i} -sub test000${i} output
output1=`avwstats output -m`
else 
avwmaths  vol00${i} -sub test00${i} output
output1=`avwstats output -m`
fi
echo $i $output1
if [ $output1 != 0.000000 ]
then
echo "Problem found in this comparision"
flag=1
fi 
i=`expr $i + 1`
done

rm vol*
rm test*
rm output*

echo ""
if [ $flag -ne 0 ]
then
echo "non-zero comparison; possible problem with avwsplit"
else 
echo "All comparisons zero. No problems detected with avwsplit"
fi