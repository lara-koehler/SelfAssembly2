# SelfAssembly2


This codes are used to run simulated annealing of particles with anisotropic interactions, coded into an interaction matrix, and analyze the results. Here is the few steps to follow before being able to use the code. 


1/ Download the files

2/ Create the lattice files
- install json and subprocess if it is not already installed in your python environemnt (with pip install for instance)
- Create a /Lattice subfolder in your working folder, where the files of the lattice geometries will be stored
- Open the file BuildLattice.py, modify the directory to the folder you just created, and run it. It should create a lot of  json files in the Lattice folder.

3/ Compile the C++ file 
Go to the terminal, go to the directory where the codes are stored and run "make". This should run the "makefile" file, it might generate a few warning but not errors. After this, most of the .hpp files now also have a .o counterpart, and an exectuable "Sys" should have been created in the same folder.
If this generates errors, it might be related to how your operating system deals with makefile, or to compilers not correctly installed on your computer. I am not a pro with that, but let me know if you cannot fix it. 

4/ Get ready to launch simulation
- Create a /Results subfolder in your working folder
- Open the file System.py, go to the class "ReadSeveralSystems" and modify the "directoryResults" to the Results folder you created.

  So far, it is the following :
  " if local :
    directoryResults = '/Users/lara/Documents/SelfAssembly2/Results/' "
  
- Open the LaunchSimulationExample jupyter notebook
- Change the "localDirectory" and "projectFolderName" to the correct ones in the first cell of the notebook. 
- Get started with the notebook!
