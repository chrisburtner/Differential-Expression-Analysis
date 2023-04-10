#! /bin/bash

# loop through sample feature tables to make individual .mat files
for table in $@
do
    # pick any one sample features table for the first step
    if [[ ! -f genes.mat ]]
    then
        echo "" > genes.mat
        cat $table | cut -f 1 >> genes.mat
    fi

    sample_name=$(basename $table | cut -f 1 -d'.')
    echo $sample_name >> $sample_name.mat
    cat $table | cut -f 2 >> $sample_name.mat
done

# paste tables together into single matrix
paste -d'\t' genes.mat [^g]*.mat > gene.counts.matrix
