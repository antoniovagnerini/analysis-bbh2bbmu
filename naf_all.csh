#!/bin/csh -f

if ( $#argv < 4 ) then
   echo Need to give sample name, the macro name, the config name and the number of splits
   exit
endif
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0- REMOVING OLD JOBS 
echo "-------------------------------------------------------------"
echo "0) WARNING: killing all OLD jobs by " $user 
qdel -u $user

set macro = $1
set config = $2
set rootfilelist = $3
set nsplit = $4

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1-SPLIT ROOTFILE LIST & SUBMIT NEW JOBS
echo "-------------------------------------------------------------"
echo "1) Splitting rootfileList & submit jobs"
./split.csh $nsplit $rootfilelist

set files = `/bin/ls *_x???.txt`

foreach file ( $files )

   set counter = `basename $file .txt | awk -F "x" '{print $2}'`
   set exedir = "NAF_"$macro"_"$counter
   if ( -d $exedir ) then
      echo "Similar jobs were already submitted. Move ot remove directories NAF_"$macro"_* and resubmit"
      exit
   endif
   mkdir -p $exedir
   cd $exedir
   mv ../$file ./rootFileList.txt
   if ( -e ../json_2017.txt ) then
      cp -p ../json_2017.txt .
   endif
   cp ../$config . 
   ../qsub.sh $macro $config
   sleep 5
   cd -
end

sleep 5m 

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2- CHECK JOB STATUS every 5min-->1h & 
echo "-------------------------------------------------------------"
echo "2) Checking jobs status up to maximum 1 h"
foreach dt ( 5m 10m 15m 20m 25m 30m 35m 40m 45m 50m 55m 60m )
    qstat > temp.txt
    if ( { grep -q $user'     r' temp.txt } ) then
	echo "After" $dt "some jobs are still Running"
	sleep 5m
	rm temp.txt
    else if ( $dt == 60m ) then
	echo "JOBs taking longer than 1 h... Aborting"
	exit 
    else
	echo "All jobs DONE after" $dt
 	break
    endif
end

rm temp.txt

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 3- OUTPUT LOGS & MERGE OUTPUT ROOTfiles
echo "-------------------------------------------------------------"
echo " 3) Saving LOG files and MERGING root output"

set NAFdirs = 'NAF_*_???'
set newdir = 'Merge_Histo'
set errdir = $newdir/ERR
set logdir = $newdir/LOG

echo "Error folder :" $errdir

rm -rf $newdir; mkdir $newdir; mkdir $errdir; mkdir $logdir

foreach NAFdir ( $NAFdirs )

    cd $NAFdir
    if ( { grep -q 'segmentation' *.out } || { grep -q 'err' *.out } || { grep -q 'nan' *.out }  ) then  # if err foud in output NAF  
	#cd  ..
	echo $NAFdir
	cat *.out | grep 'segmentation'
	head -10 *.out | grep 'error'
	cat *.out | grep 'nan'
	cp *.out ../$errdir/$NAFdir.log # copy in err folder log
	echo "-------------------------"
    else
	cp histogram*.root ../$newdir/histograms_$NAFdir.root # copy histogram_NAF_macro_???.root file
	cp *.out ../$logdir/$NAFdir.log
    endif
    cd ..; rm -rf $NAFdir
end

#Merge histograms.root files
cd $newdir
hadd histograms.root histograms_NAF_*.root
rm -rf histograms_NAF_*
exit


