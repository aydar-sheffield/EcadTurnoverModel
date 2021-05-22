This Chaste project contains source code necessary to run simulations of an extended vertex model,
where interfacial energy term depends on deviation of junctional length from rest length, and 
rest length remodels so that it tends to approach the current length. 

Additionally, the project contains the code generating multiple tissue configurations used as initial condition for the above simulations. 

To generate initial vertex meshes (i.e. tissues), run generate_meshes.sh script. This generates 12 meshes by launching 12 instances of TestGenerateRandomEpithelia. To change the number of distinct meshes, edit the script file.
The resultant meshes are saved in the test output folder used by Chaste.   

To generate simulations of the extended vertex model, run generate_many_simulations.sh script. This script launches 12 instances of TestEcadTurnover tests, each using one of the previously generated meshes as initial condition.
The test results are stored in the test output folder used by Chaste.

Prior to launching either of the above scripts, a user must configure the Chaste project as described in Chaste tutorials. 

Matlab scripts to analyse the data are also included. To run these, a user must edit the script files, specifying the paths to folders, containing simulations outputs.    
