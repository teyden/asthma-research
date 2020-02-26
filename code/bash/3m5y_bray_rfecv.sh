source ~/.bashrc

export R_LIBS=~/Projects/asthma/asthma-research/packrat/lib/x86_64-redhat-linux-gnu/3.5.1:$R_LIBS

Rscript ~/Projects/asthma/asthma-research/code/r/scripts/batch-scripts/run_rfe_cv.R --timepoint 3m --kernel bray --outcome diseasestatus_5y_binary --selectioncriterion pval
