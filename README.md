These script files use the metadata file in the compressed folder named githup_data.tar.gz deposited in the GEO repository (GSE145639) to regenerate the manuscript plots. To run 
each python or R script file, the script file takes two input arguments. The first input argument is the full path of the folder of the metadata mentioned above, for example: 
/home/user/githup_data. The second argument is the full path of the folder where to save the file of the output results, for example /home/usrer/Results.  All the script files 
mandatorily need two arguments, but some of the python script files will save the output files in the same input directory.  

Pre-requisites:
These tools/packages should be installed and defined in the system path:
Python v. 2.7
samtools
bedtools
subread
deeptools v. 3.2.1
R v. 3.5.0

