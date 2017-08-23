# ProteinStructureSearch-Simulation
##### This repository contains code for [this](https://ricsinaruto.github.io/website/docs/tdk.pdf) paper. The state of the program as described in the paper can be found in the Early-searching-algorithm branch. In order branches you can find further variants of the program and I created separate readme files for all the branches so check them out.
##### The program was basically created to easily simulate the change of dipole moments of complex protein structures and provide the user with a good looking 3D interface in OpenGL. Based on the change of dipole moments we can construct logic circuits as described in the paper mentioned before. Furthermore I implemented several searching algorithms to find protein structures that satisfy user defined logic functions.
<a><img src="https://github.com/ricsinaruto/ProteinStructureSearch-Simulation/blob/master/first_start.png" align="top" height="450" ></a>

## Get Started
##### The version of the program described and used in the [paper](https://ricsinaruto.github.io/website/docs/tdk.pdf) can be downloaded [here](https://ricsinaruto.github.io/website/docs/molecular_simulation_program.rar). You can simply extract the archive and run the executable. In order to run other versions of the program you will need to compile the project with Visual Studio. You will also need [freeglut](http://freeglut.sourceforge.net/index.php#download) library and includes.
##### You can find a detailed user guide in section 3.7 of the [paper](https://ricsinaruto.github.io/website/docs/tdk.pdf), and other instructions on how to use different parts of the program in sections 3.4.3 and 3.5.2. The user guide will mostly be true for all of the versions of the program, but if not, changes are documented in the separate branches.

## Description of Branches
The latest version of the program can be found in the GA-all-big-structures branch. The code in this version has also been extensively commented, cleaned and translated to english.
#### Early-searching-algorithm
Here you can find the original version of the program, as described in the [paper](https://ricsinaruto.github.io/website/docs/tdk.pdf). A very simple searching algorithm is used to find good protein structures. Basically it looks through all the possible protein structures and field value configurations until it finds one that satisfies the user defined logic function.
#### simulated-annealing
For this variant of the program I implemented a simulated annealing based searching algorithm to find field values to apply to a user defined protein structure in order to achieve a user defined logic function.
#### Annealing+old-algorithm
This variant of the program uses the previously mentioned simulated annealing method to search for field values, combined with the searching algorithm in the Early-searching-algorithm branch to automatically search for protein structures at the same time.
#### Annealing+random-structures
This variant of the program generates random protein structures and applies the simulated annealing algorithm to each structure to optimize the field values applied to the proteins. Similar to the previous branch, but instead of looking through all possible structures while increasing the size of the structure, we take one random sized structure at each iteration step.
#### GA-all-big-structures
In this variant of the program I implemented a genetic algorithm to find a protein structure as well as field values that satisfy a user defined logic function simultanously. This is the newest and best version of the program.
## Master folders
#### TDK-earliest version
This is kept just for archival purposes. It's the earliest instantiation of the project done with windows forms graphics.
#### TDK-generator
I used this slightly modified version of the program to quickly generate huge datasets of different protein configuration simulations. You can modify the source code in harmony_search function in interaction.cpp to generate whatever simulations you would like.
#### TDK-tester
This is a slightly modified version with which you can manually test a protein structure. Just some simple modifications to the harmony_search function in interaction.cpp.



### Credits
I used tutorial code from [here](https://github.com/davidwparker/opengl-screencasts-2) to get started with OpenGL, and to set up a basic environment.
