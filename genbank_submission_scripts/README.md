README

These scripts are to prepare files for submission to genbank.

Both take the same arguments:

--master-csv GLab master spreadsheet, downloaded in csv format. 
--file-path Path to genomic data, relative to this script 
--genbank-path Output directory where new files will be made, relative to this script
--date-submission For naming folders associated with this submission

pull_bam_files.py must be run on ruddle because that's where bam files are kept.

They pull all of the files together from individual runs and reorganise by originating lab.
We then submit by originating lab.
