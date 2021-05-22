#!/bin/bash
#Generates N initial meshes
counter=11
N=12
while [ $counter -le $N ]
do
	../../../build/projects/CompToolsExamples/test/TestGenerateRandomEpithelia -opt $counter
	echo $counter
	((counter++))
done
echo All done
