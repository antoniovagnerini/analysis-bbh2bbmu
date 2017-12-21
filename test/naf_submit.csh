#!/bin/csh -f

if ( $#argv < 4 ) then
   echo Need to give sample name, the macro name, the config name and the number of splits
   exit
endif

set macro = $1
set config = $2
set rootfilelist = $3
set nsplit = $4

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


exit




#./qsub.sh $exeName $sampleName
