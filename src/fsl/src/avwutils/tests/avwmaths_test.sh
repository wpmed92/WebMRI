#!/bin/sh
#testing avwmaths
flag=0
#avwmaths volume test
for setting in -T -X -Y -Z 
do
for loop in mean std max maxn median
do
set2=`echo $setting$loop`
../avwmaths_32R $FSLTESTDIR/common/filtered_func_data.nii.gz $set2 temp
../avwmaths++ $FSLTESTDIR/common/filtered_func_data.nii.gz $set2 temp1 -odt float
avwmaths temp1 -Tmean temp2 
avwmaths temp -Tmean temp1
avwmaths temp1 -sub temp2 -abs temp
output1=`avwstats temp -m`
echo "Result of $setting$loop test:" $output1
if [ $output1 != 0.000000 ];
then
flag=1
echo "Possible problem with $setting$loop test in avwmaths++";
fi
rm temp*
done
set2=`echo ${setting}perc`
avwmaths_32R $FSLTESTDIR/common/filtered_func_data.nii.gz $set2 40 temp
../avwmaths++ $FSLTESTDIR/common/filtered_func_data.nii.gz $set2 40 temp1 -odt float
avwmaths temp1 -Tmean temp2 
avwmaths temp -Tmean temp1
avwmaths temp1 -sub temp2 -abs temp
output1=`avwstats temp -m`
echo "Result of $set2 test:" $output1
if [ $output1 != 0.000000 ];
then
flag=1
echo "Possible problem with $set2 test in avwmaths++";
fi
rm temp*
set2=`echo ${setting}ar1`
set3=`echo ${setting}mean`
../avwmaths++ $FSLTESTDIR/common/filtered_func_data.nii.gz $set3 tempmean -odt float
if [ $set2 = -Xar1 ]
then
avwmerge -x tempmean tempmean tempmean #2
avwmerge -x tempmean tempmean tempmean #4
avwmerge -x tempmean tempmean tempmean #8
avwmerge -x tempmean tempmean tempmean #16
avwmerge -x tempmean tempmean tempmean tempmean tempmean #32
fi
if [ $set2 = -Yar1 ]
then
avwmerge -y tempmean tempmean tempmean #2
avwmerge -y tempmean tempmean tempmean #4
avwmerge -y tempmean tempmean tempmean #8
avwmerge -y tempmean tempmean tempmean #16
avwmerge -y tempmean tempmean tempmean tempmean tempmean #32
fi
if [ $set2 = -Zar1 ]
then
../avwmerge++ -z tempout  tempmean tempmean tempmean #3
../avwmerge++ -z tempmean tempout tempout tempout tempout tempout tempout tempout #21
fi
../avwmaths++ $FSLTESTDIR/common/filtered_func_data.nii.gz -sub tempmean tempmean -odt float
../avwmaths_32R tempmean $set2 temp 
../avwmaths++  $FSLTESTDIR/common/filtered_func_data.nii.gz $set2 temp1 -odt float

../avwmaths++ temp1 -Tmean temp2 
avwmaths temp -Tmean temp1
avwmaths temp1 -sub temp2 -abs temp
output1=`avwstats temp -m`
echo "Result of $set2 test:" $output1
if [ $output1 != 0.000000 ];
then
flag=1
echo "Possible problem with $set2 test in avwmaths++";
fi
rm temp*

done













ip_32R $FSLTESTDIR/common/filtered_func_data.nii.gz temp 0 -t 10 60
../avwmaths++ $FSLTESTDIR/common/filtered_func_data.nii.gz -bptf 10 60 temp1 -odt float
avwmaths temp1 -Tmean temp2 
avwmaths temp -Tmean temp1
avwmaths temp1 -sub temp2 -abs temp
output1=`avwstats temp -m`
echo "Result of bptf test:" $output1
limit=0.05  #Good value, test seems sensitive to changes in algorithm/input errors
output2=`echo $output1 $limit | awk '{ print ($1 < $2) ? "passed" : "failed" }' `
echo "$output2"
if [ $output2 != passed ];
then
flag=1
echo "Possible problem with bptf test in avwmaths++";
else echo "bptf test within tolerance";
fi
rm temp*

avwmaths $FSLTESTDIR/common/filtered_func_data.nii.gz -Tmean tempvol
../avwconv  --dilate -b 10 -i tempvol -o temp
../avwmaths++ tempvol -kernel box 10 -dilF temp1 
avwmaths temp1 -Tmean temp2 
avwmaths temp -Tmean temp1
avwmaths temp1 -sub temp2 -abs temp
output1=`avwstats temp -m`
echo "Result of avwconv erosion test:" $output1
limit=0.0001
output2=`echo $output1 $limit | awk '{ print ($1 < $2) ? "passed" : "failed" }' `
echo "$output2"
if [ $output2 != passed ];
then
flag=1
echo "Possible problem with avwconv erosion test in avwmaths++";
else echo "avwconv erosion test within tolerance";
fi
rm temp*

avwmaths $FSLTESTDIR/common/filtered_func_data.nii.gz -Tmean tempvol
../avwconv  --erode -b 10 -i tempvol -o temp
../avwmaths++ tempvol -kernel box 10 -eroF temp1 
avwmaths temp1 -Tmean temp2 
avwmaths temp -Tmean temp1
avwmaths temp1 -sub temp2 -abs temp
output1=`avwstats temp -m`
echo "Result of avwconv dilation edge test:" $output1
limit=0.0001
output2=`echo $output1 $limit | awk '{ print ($1 < $2) ? "passed" : "failed" }' `
echo "$output2"
if [ $output2 != passed ];
then
flag=1
echo "Possible problem with avwconv dilation test in avwmaths++";
else echo "avwconv dilation test within tolerance";
fi
rm temp*

avwmaths $FSLTESTDIR/common/filtered_func_data.nii.gz -Tmean tempvol
avwconv  -b 40 -i tempvol -o temp
../avwmaths++ tempvol -kernel box 40 -fmeanu temp1 -odt float
avwmaths temp1 -Tmean temp2 
avwmaths temp -Tmean temp1
avwmaths temp1 -sub temp2 -abs temp
output1=`avwstats temp -m`
echo "Result of large box convolve edge test:" $output1
limit=0.0001
output2=`echo $output1 $limit | awk '{ print ($1 < $2) ? "passed" : "failed" }' `
echo "$output2"
if [ $output2 != passed ];
then
flag=1
echo "Possible problem with large box convolve test in avwmaths++";
else echo "large box convolve test within tolerance";
fi
rm temp*

avwmaths $FSLTESTDIR/common/filtered_func_data.nii.gz -Tmean tempvol
../avwconv  -b 8 -i tempvol -o temp
../avwmaths++ tempvol -kernel box 8 -fmeanu temp1 -odt float
avwmaths temp1 -Tmean temp2 
avwmaths temp -Tmean temp1
avwmaths temp1 -sub temp2 -abs temp
output1=`avwstats temp -m`
echo "Result of small box convolve edge test:" $output1
limit=0.0001
output2=`echo $output1 $limit | awk '{ print ($1 < $2) ? "passed" : "failed" }' `
echo "$output2"
if [ $output2 != passed ];
then
flag=1
echo "Possible problem with small box convolve test in avwmaths++";
else echo "small box convolve test within tolerance";
fi
rm temp*

avwmaths $FSLTESTDIR/common/filtered_func_data.nii.gz -Tmean tempvol
../avwconv  -s 12 --median -i tempvol -o temp
../avwmaths++ tempvol -kernel sphere 12 -fmedian temp1 
avwmaths temp1 -Tmean temp2 
avwmaths temp -Tmean temp1
avwmaths temp1 -sub temp2 -abs temp
output1=`avwstats temp -m`
echo "Result of spherical kernel test:" $output1
limit=0.0001
output2=`echo $output1 $limit | awk '{ print ($1 < $2) ? "passed" : "failed" }' `
echo "$output2"
if [ $output2 != passed ];
then
flag=1
echo "Possible problem with spherical kernel test in avwmaths++";
else echo "spherical kernel test within tolerance";
fi
rm temp*

ip_16SI $FSLTESTDIR/common/filtered_func_data.nii.gz temp 0 -S 3
../avwmaths++ $FSLTESTDIR/common/filtered_func_data.nii.gz -kernel box 12 -fmedian temp1
avwmaths temp1 -Tmean temp2 
avwmaths temp -Tmean temp1
avwmaths temp1 -sub temp2 -abs temp
output1=`avwstats temp -m`
echo "Result of ip median box kernel test:" $output1
if [ $output1 != 0.000000 ];
then
flag=1
echo "Possible problem with  ip median box kernel test in avwmaths++";
fi
rm temp*

ip_32R $FSLTESTDIR/common/filtered_func_data.nii.gz temp 0 -s 4
../avwmaths++ $FSLTESTDIR/common/filtered_func_data.nii.gz -kernel gauss 4 -fmean temp1 -odt float
avwmaths temp1 -Tmean temp2 
avwmaths temp -Tmean temp1
avwmaths temp1 -sub temp2 -abs temp
output1=`avwstats temp -m`
echo "Result of ip gaussian kernel test:" $output1
limit=0.0001
output2=`echo $output1 $limit | awk '{ print ($1 < $2) ? "passed" : "failed" }' `
echo "$output2"
if [ $output2 != passed ];
then
flag=1
echo "Possible problem with ip gaussian kernel test in avwmaths++";
else echo "ip gaussian kernel test within tolerance";
fi
rm temp*

avwmaths_32R $FSLTESTDIR/common/filtered_func_data.nii.gz -roi 10 30 10 30 10 30 0 10 temp
../avwmaths++ $FSLTESTDIR/common/filtered_func_data.nii.gz -roi 10 30 10 30 10 30 0 10 temp1 -odt float
avwmaths temp1 -Tmean temp2 
avwmaths temp -Tmean temp1
avwmaths temp1 -sub temp2 -abs temp
output1=`avwstats temp -m`
echo "Result of roi test:" $output1
if [ $output1 != 0.000000 ];
then
flag=1
echo "Possible problem with roi test in avwmaths++";
fi
rm temp*

ip_32R $FSLTESTDIR/common/filtered_func_data.nii.gz temp 90 -m mask1
../avwmaths++  $FSLTESTDIR/common/filtered_func_data.nii.gz -thrp 90 -Tmin -bin mask -odt float
../avwmaths++ $FSLTESTDIR/common/filtered_func_data.nii.gz -mas mask temp1 -odt float
avwmaths temp1 -Tmean temp2 
avwmaths temp -Tmean temp1
avwmaths temp1 -sub temp2 -abs temp
output1=`avwstats temp -m`
echo "Result of ip thresholding test:" $output1
limit=1
output2=`echo $output1 $limit | awk '{ print ($1 < $2) ? "passed" : "failed" }' `
echo "$output2"
if [ $output2 != passed ];
then
flag=1
echo "Possible problem with ip thresholding (-thrmt) test in avwmaths++";
else echo "ip thresholding test within tolerance";
fi
rm temp*
rm mask*


for loop in -sub -add -mul -div -mas -max -min
do
avwmaths_32R $FSLTESTDIR/common/filtered_func_data.nii.gz $loop $FSLTESTDIR/common/mask.nii.gz temp
../avwmaths++ $FSLTESTDIR/common/filtered_func_data.nii.gz $loop $FSLTESTDIR/common/mask.nii.gz temp1 -odt float
avwmaths temp1 -Tmean temp2 
avwmaths temp -Tmean temp1
avwmaths temp1 -sub temp2 -abs temp
output1=`avwstats temp -m`
echo "Result of $loop volume test:" $output1
if [ $output1 != 0.000000 ];
then
flag=1
echo "Possible problem with $loop volume test in avwmaths++";
fi
rm temp*
done


for loop in -exp -log -sqr -sqrt -abs -edge -bin
do
avwmaths_32R $FSLTESTDIR/common/filtered_func_data.nii.gz $loop temp
../avwmaths++ $FSLTESTDIR/common/filtered_func_data.nii.gz $loop temp1 -odt float
avwmaths temp1 -Tmean temp2 
avwmaths temp -Tmean temp1
avwmaths temp1 -sub temp2 -abs temp
output1=`avwstats temp -m`
echo "Result of $loop volume test:" $output1
if [ $output1 != 0.000000 ];
then
flag=1
echo "Possible problem with $loop volume test in avwmaths++";
fi
rm temp*
done

../avwmaths++  $FSLTESTDIR/common/filtered_func_data.nii.gz -exp infvol -odt float
for loop in -nan -nanm 
do
avwmaths_32R infvol $loop temp
../avwmaths++ infvol $loop temp1 -odt float
avwmaths temp1 -Tmean temp2 
avwmaths temp -Tmean temp1
avwmaths temp1 -sub temp2 -abs temp
output1=`avwstats temp -m`
echo "Result of $loop volume test:" $output1
if [ $output1 != 0.000000 ];
then
flag=1
echo "Possible problem with $loop volume test in avwmaths++";
fi
rm temp*
done
rm infvol*

ip_32R $FSLTESTDIR/common/filtered_func_data.nii.gz temp 0 -i 2000
../avwmaths++ $FSLTESTDIR/common/filtered_func_data.nii.gz -inm 2000 temp1 -odt float
avwmaths temp1 -Tmean temp2 
avwmaths temp -Tmean temp1
avwmaths temp1 -sub temp2 -abs temp
output1=`avwstats temp -m`
echo "Result of ip norm test:" $output1
limit=4
output2=`echo $output1 $limit | awk '{ print ($1 < $2) ? "passed" : "failed" }' `
echo "$output2"
if [ $output2 != passed ];
then
flag=1
echo "Possible problem with ip norm (-inm) test in avwmaths++";
else echo "ip norm test within tolerance";
fi
rm temp*

ip_32R $FSLTESTDIR/common/filtered_func_data.nii.gz temp 0 -I 2000
../avwmaths++ $FSLTESTDIR/common/filtered_func_data.nii.gz -ing 2000 temp1 -odt float
avwmaths temp1 -Tmean temp2 
avwmaths temp -Tmean temp1
avwmaths temp1 -sub temp2 -abs temp
output1=`avwstats temp -m`
echo "Result of ip global norm test:" $output1
limit=4
output2=`echo $output1 $limit | awk '{ print ($1 < $2) ? "passed" : "failed" }' `
echo "$output2"
if [ $output2 != passed ];
then
flag=1
echo "Possible problem with ip global norm (-ing) test in avwmaths++";
else echo "ip global norm test within tolerance";
fi
rm temp*

#avwmaths data tests
for loop in -sub -add -mul -div -thr -uthr -max -min
do
avwmaths_32R $FSLTESTDIR/common/filtered_func_data.nii.gz $loop 10 temp
../avwmaths++ $FSLTESTDIR/common/filtered_func_data.nii.gz $loop 10 temp1 -odt float
avwmaths temp1 -Tmean temp2 
avwmaths temp -Tmean temp1
avwmaths temp1 -sub temp2 -abs temp
output1=`avwstats temp -m`
echo "Result of $loop data test:" $output1
if [ $output1 != 0.000000 ];
then
flag=1
echo "Possible problem with $loop data test in avwmaths++";
fi
rm temp*
done
rm kernel*




echo ""
if [ $flag -ne 0 ]
then
echo "non-zero comparison; possible problem with avwmaths++"
else 
echo "All comparisons zero. No problems detected with avwmaths++"
fi