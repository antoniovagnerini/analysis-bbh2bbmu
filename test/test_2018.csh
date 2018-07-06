#!/bin/csh -f

set datastream = "PROMPTRECO"
set macro = "SimpleMssmHbbmuAnalysis"
set template = "template_csv"                                                                                                                                                                          


set outdir = "ROOTFILES/"$datastream
if (! -d $outdir ) then
    echo Output dir $outdir created
    mkdir $outdir
endif

cp ROOTFILELIST/$datastream/json_2018.txt json_2018.txt
                                                                                                                              
foreach era ( 2018A )                                                                                                                                                                              
    echo Processing ERA $era
    cp ROOTFILELIST/$datastream/rootFileList_$era.txt rootFileList.txt 
    ./htc_all.csh $macro $template.cfg rootFileList.txt 5
    mv Merge_Histo/histograms.root $outdir/histograms_${era}csv.root
    echo '================================================'
end

set template = "template_deepcsv" 

foreach era ( 2018A )
    echo Processing ERA $era
    cp ROOTFILELIST/$datastream/rootFileList_$era.txt rootFileList.txt
    ./htc_all.csh $macro $template.cfg rootFileList.txt 5
    mv Merge_Histo/histograms.root $outdir/histograms_${era}deepcsv.root                                                                                                                                 
    echo '================================================'
end

'''
