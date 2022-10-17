#!/bin/sh
#test 3D volume
flag=0 
for setting in -t -x -y -z 
do
avwmerge $setting temp $FSLTESTDIR/common/filtered_func_data.nii.gz $FSLTESTDIR/common/filtered_func_data2.nii.gz 
../avwmerge++ $setting temp2 $FSLTESTDIR/common/filtered_func_data.nii.gz $FSLTESTDIR/common/filtered_func_data2.nii.gz 
avwmaths temp2 -Tmean temp3 
avwmaths temp -Tmean temp2
avwmaths temp3 -sub temp2 temp
output1=`avwstats temp -m`
echo "Result of $setting merge test:" $output1
if [ $output1 != 0.000000 ];
then
flag=1
echo "Possible problem with $setting test in avwmerge++";
fi
rm temp*
done

avwmerge -a temp $FSLTESTDIR/common//vol*
../avwmerge++ -a temp2 $FSLTESTDIR/common/vol*
avwmaths temp2 -Tmean temp3 
avwmaths temp -Tmean temp2
avwmaths temp3 -sub temp2 temp
output1=`avwstats temp -m`

echo "Result of -a merge test:" $output1
if [ $output1 != 0.000000 ];
then
flag=1
echo "Possible problem with -a test in avwmerge++";
fi

rm temp*






echo ""
if [ $flag -ne 0 ]
then
echo "non-zero comparison; possible problem with avwmerge++"
else 
echo "All comparisons zero. No problems detected with avwmerge++"
fi