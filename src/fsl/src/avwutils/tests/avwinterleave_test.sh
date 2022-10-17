#!/bin/sh
flag_maths=0 
flag_bit=0
avwinterleave $FSLTESTDIR/common/vol0000.nii.gz  $FSLTESTDIR/common/vol0001.nii.gz temp1 
../avwinterleave++  $FSLTESTDIR/common/vol0000.nii.gz  $FSLTESTDIR/common/vol0001.nii.gz temp2 
cmp temp1.nii.gz temp2.nii.gz
output2=$?
if [ $output2 -ne 0 ]
then
flag_bit=1
echo "Standard merge failed bit-comparison test in avwinterleave++"
fi
avwmaths temp1 -sub temp2 temp3
output1=`avwstats temp3 -m`
echo "Result of interleave test:" $output1  $output2
if [ $output1 != 0.000000 ];
then
flag_maths=1
echo "standard merge failed maths test in avwinterleave++";
fi

rm temp*

avwinterleave $FSLTESTDIR/common/vol0000.nii.gz  $FSLTESTDIR/common/vol0001.nii.gz temp1 -i
../avwinterleave++  $FSLTESTDIR/common/vol0000.nii.gz  $FSLTESTDIR/common/vol0001.nii.gz temp2 -i
cmp temp1.nii.gz temp2.nii.gz
output2=$?
if [ $output2 -ne 0 ]
then
flag_bit=1
echo "inverse merge failed bit-comparison test in avwinterleave++"
fi
avwmaths temp1 -sub temp2 temp3
output1=`avwstats temp3 -m`
echo "Result of inverse interleave test:" $output1  $output2
if [ $output1 != 0.000000 ];
then
flag_maths=1
echo "inverse merge failed maths test in avwinterleave++";
fi

rm temp*

echo ""
if [ $flag_maths -ne 0 -o $flag_bit -ne 0 ]
then
echo "non-zero comparison; possible problem with avwinterleave++"
else 
echo "All comparisons zero. No problems detected with avwinterleave++"
fi