#!/bin/bash

counter=0
N=11
echo $counter
for ((i=0;i<=N;i+=1)); do
	../../../build/projects/CompToolsExamples/test/TestEcadTurnover -opt $i &
done
echo All done
