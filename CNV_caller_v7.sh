#!/bin/bash

# Wrapper for submitting CNV to Medivation cluster. This script
# will create a .pbs script for each individual job requested.

export OUT="/work/trqlogs/cfDNA/"

samplelist=$1

while IFS=$'\t' read -r line; do
    bam=`echo $line | awk '{print$1}'`
    name=`echo $line | awk '{print$2}'`
    output=`echo $line | awk '{print$3}'`
    background=`echo $line | awk '{print $4}'`
    runexac=`echo $line | awk '{print $5}'`
    echo $name
    echo $bam
    echo $output 
    cat > $OUT$name".pbs" << EOF

#!/bin/bash
#
#This is submission script for torque queue
# it runs CNV.R code each each sample and .bam
# file should be placed in SUBFILE folder
# 
#These commands set up the Grid Environment for your job:
#PBS -j oe
#PBS -N $name
#PBS -l nodes=1:ppn=26,walltime=10:00:00
#PBS -V
#PBS -o localhost:/work/trqlogs/cfDNA/out.$name
#PBS -e localhost:/work/trqlogs/cfDNA/error.$name
#cPBS -m a 
#cPBS -koe

SRCDIR=/tools/cfDNA/bin
RAWDIR=$bam
SUBFILE="$name";
WORKDIR=$output
BCK=$background
RE=$runexac
EOF

    cat >> $OUT$name".pbs" <<'EOF'

### actual job: (example: Rscript CNV_v4.R /path/to/bam/test.bam samplename /path/to/output)
Rscript $SRCDIR/CNV_v7.R $RAWDIR $SUBFILE $WORKDIR $BCK $RE

sleep 2
EOF

eval `qsub $OUT$name".pbs"`
done < "$samplelist"    
