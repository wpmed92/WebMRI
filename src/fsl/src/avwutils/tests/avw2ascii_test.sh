#!/bin/sh
#test 3D volume
avw2ascii $FSLTESTDIR/common/filtered_func_data.nii.gz test
../avw2ascii++ $FSLTESTDIR/common/filtered_func_data.nii.gz test2
flag=0
i=0
while [ $i -lt 100 ] 
do
if [ $i -lt 10 ] 
then 
cmp test0000${i} test20000${i} 
else 
cmp test000${i} test2000${i}
fi
j=$?
echo $i $j
if [ $j -ne 0 ]
then
echo "Problem found in this comparision"
flag=1
fi 
i=`expr $i + 1`
done

rm test*

echo ""
if [ $flag -ne 0 ]
then
echo "non-zero comparison; possible problem with avw2ascii++"
else 
echo "All comparisons zero. No problems detected with avw2ascii++"
fi