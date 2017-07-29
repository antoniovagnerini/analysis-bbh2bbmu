#!/bin/csh -f

if ( $#argv < 2 ) then
   echo "Usage: split.csh num_lines file"
   exit
endif

set nlines = $1 
set infile = $2

set filename = `basename $infile | awk -F "." '{print $1}'`
set suffix   = `basename $infile | awk -F "." '{print $2}'`

split -a 3 -d -l $nlines $infile

set files = `/bin/ls -1 x???`

foreach f ( $files )
   mv $f $filename"_"$f".txt"
end

exit
