#!/bin/sh 
# $1 - code name
cat > $1.zsh <<EOF
#!/bin/zsh
#
#(make sure the right shell will be used)
#$ -S /bin/zsh
#
#(the cpu time for this job)
#$ -l h_cpu=7:59:00
#
#(the maximum memory usage of this job)
#$ -l h_vmem=2G
#
#(use hh site)
#$ -l site=hh 
#(stderr and stdout are merged together to stdout)
#$ -j y
#
# use SL6
#$ -l os=sld6
#
# use current dir and current environment
#$ -cwd
#$ -V
#
#$ -o $1.out
#
#$ -e $1.err
$1

EOF

chmod u+x $1.zsh
qsub -q "short.q" $1.zsh
