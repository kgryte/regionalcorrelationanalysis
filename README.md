Regional Correlation Analysis
===========================

This program implements 2D Regional Correlation Analysis, an algorithm proposed in the paper:


Xuefei Wang, H. Peter Lu (2008) "2D Regional Correlation Analysis of Single-Molecule Trajectories," J. Phys. Chem. B. 112, 1492-14926



## Compilation Notes


* Need to compile using -lm option so sqrt() can take variable input.

* LINUX:
	* ...to create the executable...
  	 	- g++ -o regionalCorrelationAnalysis regionalCorrelationAnalysis.cpp -lm
	* ...to run...
 		- ./regionalCorrelationAnalysis <inputDataFilename> <outputFilename>
 			e.g., ./regionalCorrelationAnalysis intputData.dat output.dat
 


## History


2012-06-25: KGryte. Created.



## Contact


Contact:

Kristofer Gryte  
University of Oxford  
Oxford, UK  
k.gryte1(at)physics.ox.ac.uk
  
Copyright (c) 2012-2013







