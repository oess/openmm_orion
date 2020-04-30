This Floe is for Gromacs MD experts who simply want to run Gromacs using as
input their own (uploaded to Orion) Gromacs *.tpr* files. The floe will run 
for *n* hours (10hrs default), outputting the
trajectory file and a recovery dataset. The Gromacs cube 
will then restart from the recovery dataset and run for
an additional *n* hours in a cycle
until the number of md steps specified in 
the *.tpr* file are finished. If the recovery dataset is
provided as input, Gromacs will start from the last
checkpoint saved in the recovery dataset. If both *.tpr* 
and recovery dataset files are provided the recovery dataset 
will overwrite the *.tpr* file.
