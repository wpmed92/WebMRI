#!/bin/sh
flag=0 
cp  $FSLTESTDIR/common/filtered_func_data.nii.gz temp.nii.gz
cp  $FSLTESTDIR/common/filtered_func_data.nii.gz temp2.nii.gz
avwcpgeom $FSLTESTDIR/common/avg152T1_brain.hdr temp.nii.gz
../avwcpgeom++ $FSLTESTDIR/common/avg152T1_brain.hdr temp2.nii.gz

cmp temp.nii.gz temp2.nii.gz
j=$?
echo $j
if [ $j -ne 0 ]
then
echo "Problem found in 4D:3D comparison in avwcpgeom++"
flag=1
fi 

rm temp*

cp  $FSLTESTDIR/common/filtered_func_data.nii.gz temp.nii.gz
cp  $FSLTESTDIR/common/filtered_func_data.nii.gz temp2.nii.gz
avwcpgeom $FSLTESTDIR/common/avg152T1_brain.hdr temp.nii.gz -d
../avwcpgeom++ $FSLTESTDIR/common/avg152T1_brain.hdr temp2.nii.gz -d

cmp temp.nii.gz temp2.nii.gz
j=$?
echo $j
if [ $j -ne 0 ]
then
echo "Problem found in 4D:3D comparison in avwcpgeom++"
flag=1
fi 

rm temp*

cp  $FSLTESTDIR/common/avg152T1_brain.hdr temp.hdr
cp  $FSLTESTDIR/common/avg152T1_brain.img temp.img
cp  $FSLTESTDIR/common/avg152T1_brain.hdr temp2.hdr
cp  $FSLTESTDIR/common/avg152T1_brain.img temp2.img

avwcpgeom  $FSLTESTDIR/common/filtered_func_data.nii.gz temp
../avwcpgeom++  $FSLTESTDIR/common/filtered_func_data.nii.gz temp2

cmp temp.hdr temp2.hdr
j=$?
echo $j
if [ $j -ne 0 ]
then
echo "Problem found in 3D:4D comparison in avwcpgeom++"
flag=1
fi 

cmp temp.img temp2.img
j=$?
echo $j
if [ $j -ne 0 ]
then
echo "Problem found in 3D:4D comparison in avwcpgeom++"
flag=1
fi 

rm temp*

cp  $FSLTESTDIR/common/avg152T1_brain.hdr temp.hdr
cp  $FSLTESTDIR/common/avg152T1_brain.img temp.img
cp  $FSLTESTDIR/common/avg152T1_brain.hdr temp2.hdr
cp  $FSLTESTDIR/common/avg152T1_brain.img temp2.img

avwcpgeom  $FSLTESTDIR/common/filtered_func_data.nii.gz temp -d
../avwcpgeom++  $FSLTESTDIR/common/filtered_func_data.nii.gz temp2 -d

cmp temp.hdr temp2.hdr
j=$?
echo $j
if [ $j -ne 0 ]
then
echo "Problem found in 3D:4D comparison in avwcpgeom++"
flag=1
fi 

cmp temp.img temp2.img
j=$?
echo $j
if [ $j -ne 0 ]
then
echo "Problem found in 3D:4D comparison in avwcpgeom++"
flag=1
fi 

rm temp*

echo ""
if [ $flag -ne 0 ]
then
echo "non-zero comparison; possible problem with avwcpgeom++"
else 
echo "All comparisons zero. No problems detected with avwcpgeom++"
fi