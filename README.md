# Hug
Hug model : a tool to analyse geological mixing system

# Introduction

This file contains explanations concerning the data and the programs used to obtain the results presented in the paper "An interaction point process for Bayesian inference of the composition of multiple water source from hydrochemical data" by C. Reype, R.S. Stoica, A. Richard, M. Deaconu.

The software and data files in this package are free for use and modification for public research and public academic activities only. This software is continously under development, and although we have done our best to ensure that the programs work correctly, we can not guarantee that they are bug free. So the responsibility for the results and conclusions obtained by using thissoftware falls to the user.

For more information feel free to contact Christophe Reype at: chris.reype@gmail.com .

The archive contains the following directories:

C++/SRC : C++ source code of the programs.
C++/DATA : data files used in the paper (synthetic, Mexico, Athabasca datasets).
C++/EXEC : contains the batch files that demonstrate how to run the programs. 
C++/PARAMS : parameter files for the executable programs.
C++/RESULTS : results files used in the paper.
R : R code for the preparation and the presentation of the data. 

# 1. DATA
The data files used can be found in the directory C++/DATA. This directory contains the synthetic data set,the Mexico data set and the Athabasca data set.

# 2. PROGRAMS
	The program sources are in the directory C++/SRC. The programs that we have written to implement our model are built using the standard C++ object oriented programming language. The programs are derived from DRLib, two C++ library. The DRLib is used to modeling, simulation and inference of markedpoint processes and is developped at Institut Elie Cartan de Lorraine by D. Gemmerlé, R.S. Stoica, C. Reype and N. Gillot. The associated documentation is available at https://gitlab.univ-lorraine.fr/labos/iecl/drlib.
	
	We have written no extra documentation For the present programs. But the programs themselves are written in explicit style, giving clear names to the variables and including many comments. The user can compile the programs using the provided Makefile. 

	The directory C++/EXEC contains the batch files showing how the C++ programs have to be used. 
	
	The directory C++/PARAMS contains the parameter files needed by the programs, and referred to in the batch files of the EXEC directory.
	
	The directory R contains a script, "plot_Hug.R", which can be used to plot the results. It also includes an executable, "Hug_application.bat",which can be used to normalise the data, create the region of interest, creates the parameters file and plotting the results. To work properly, the folder "Rshiny.7z" must be extracted into this directory.
	
	The user has to check that the path to the data file is correct.
