#!/bin/csh -f

if ( $#argv < 3 ) then
   echo Need to give sample name, the macro name and the number of splits
   exit
endif

set rootfilelist = $1
set macro = $2
set nsplit = $3


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
   if ( -e ../json.txt ) then
      cp -p ../json.txt .
   endif
   ../qsub.sh $macro
   sleep 5
   cd -
end


exit




#./qsub.sh $exeName $sampleName
