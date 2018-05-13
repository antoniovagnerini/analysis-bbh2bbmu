#!/bin/csh -f

set datastream = "BTAGCSV/RERECO"
set macro = "SimpleMssmHbbmuAnalysis"
set template = "template"

set outdir = "ROOTFILES/"$datastream

if (! -d $outdir ) then
    echo Output dir $outdir created
    mkdir $outdir
endif

cp ROOTFILELIST/$datastream/json_2017.txt json_2017.txt

foreach era ( C D E CtoE F )                                                                                                                                
    echo Processing ERA 2017$era
    cp ROOTFILELIST/$datastream/rootFileList_2017$era.txt rootFileList.txt   
    ./naf_all.csh $macro $template.cfg rootFileList.txt 5
    mv Merge_Histo/histograms.root $outdir/histograms_2017$era.root
    echo '================================================'
end


