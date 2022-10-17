#!/bin/sh

echo "************************************"
echo "Running avwroi tests"
./avwroi_test.sh
echo "************************************"
echo "Running avwnvols tests"
./avwnvols_test.sh
echo "************************************"
echo "Running avw2ascii tests"
./avw2ascii_test.sh
echo "************************************"
echo "Running avwsplit tests"
./avwsplit_test.sh
echo "************************************"
echo "Running avwcc tests"
./avwcc_test.sh
echo "************************************"
echo "Running avwmerge tests"
./avwmerge_test.sh
echo "************************************"
echo "Running avwinterleave tests"
./avwinterleave_test.sh
echo "************************************"
echo "Running avwcpgeom tests"
./avwcpgeom_test.sh
echo "************************************"
echo "Running avwhd tests"
./avwhd_test.sh
echo "************************************"
echo "Running avwcreatehd tests"
./avwcreatehd_test.sh
echo "************************************"
echo "Running avwmaths tests"
./avwmaths_test.sh
echo "************************************"