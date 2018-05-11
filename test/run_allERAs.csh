#!/bin/csh -f

set datastream = "RERECO"
#set datastream = "PROMPTRECO"
set outdir = "ROOTFILES/"$datastream

echo $outdir

if (! -d $outdir ) then
    echo Output dir $outdir created
    mkdir $outdir
endif

cp ROOTFILELIST/$datastream/json_2017.txt json_2017.txt

foreach era ( C D E F CtoE all )                                                                                                                                
    echo Processing ERA 2017$era
    cp ROOTFILELIST/$datastream/rootFileList_2017$era.txt rootFileList.txt   
    ./naf_all.csh SimpleMssmHbbmuAnalysis template.cfg rootFileList.txt 5
    mv Merge_Histo/histograms.root $outdir/histograms_2017$era.root
    echo '================================================'
end

