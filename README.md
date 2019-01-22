# ccs-post-filters

Original filters by Allison, I modified the script into 3 sub-scripts to run it in an automatic way for the Exome-Recap project.


Runs extra filters on Roslin 4.1 analysis mafs, to apply strand bias and mapping quality filters 


How to run:

Rscript <....Step1.R> #fillouts -- needs analysis mafs + bams + sample-pairing

Rscript <....Step2.R> #vcfs -- needs analysis mafs + vcfs

Rscript <....Step3.R> #apply filters -- needs analysis mafs
