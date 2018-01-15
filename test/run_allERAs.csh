#!/bin/csh -f

set outdir = "ROOTFILES"

if (! -d $outdir ) then
    echo Output dir $outdir created
    mkdir $outdir
endif

foreach era ( C-v1 C-v2 C-v3 D E )
    echo Processing ERA 2017$era
    cp ROOTFILELIST/rootFileList_2017$era.txt rootFileList.txt
    ./naf_all.csh SimpleMssmHbbmuAnalysis template.cfg rootFileList.txt 3
    mv Merge_Histo/histograms.root $outdir/histograms_2017$era.root
    #mv Merge_Histo STORE/2017$era
    echo '================================================'
end

