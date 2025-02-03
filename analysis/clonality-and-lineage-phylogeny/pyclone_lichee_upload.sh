#!/bin/bash
pyclone-vi write-results-file -i /pathtofile/results_dir/hd5/pycloneOutputFile.hd5 -o /pathtofile/results_dir/tsv/pycloneOutputFile.tsv
java -jar LICHeE/release/lichee.jar -build -i /pathtofile/results_dir/lichee_inputs/lichee_sSNV.txt -o /pathtofile/results_dir/lichee_outputs/lichee.trees.txt -cp -sampleProfile -n 4 -clustersFile /pathtofile/results_dir/lichee_inputs/lichee_cluster.txt -color -dot
