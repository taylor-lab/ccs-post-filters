# ccs-post-filters

Original filters by Allison, I modified the script into 3 sub-scripts to run it in an automatic way for the Exome-Recap project.


Runs extra filters on Roslin 4.1 analysis mafs, to apply strand bias and mapping quality filters 


How to run:

Rscript <....Step1.R> #fillouts -- needs bams and analysis mafs

Rscript <....Step2.R> #vcfs -- needs sample mapping and vcfs

Rscript <....Step3.R> #apply filters
