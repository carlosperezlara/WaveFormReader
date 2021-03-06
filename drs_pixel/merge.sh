!/bin/sh
run=${1}
echo $run
fileDRS="../physics_data/Run_${run}.dat"
filePix="../pixel_data/Run${run}_CMSTiming_FastTriggerStream_converted.root"
fileOut="Run_${run}.root"
echo $fileDRS
echo $filePix
echo $fileOut
cp $fileDRS ./drs.dat
cp $filePix ./pix.root
./maketree --inputFileName=drs.dat --pixelInputFileName=pix.root --outputFileName=out.root --nEvents=1000000
rm drs.dat pix.root
mv ./out.root merged/$fileOut
