# ProteinStructureSearch-Simulation
##### This repository contains code for [this](https://ricsinaruto.github.io/website/docs/tdk.pdf) paper. The state of the program as described in the paper can be found in the Early-searching-algorithm branch. In order branches you can find further variants of the program and I created separate readme files for all the branches so check them out.
##### The program was basically created to easily simulate the change of dipole moments of complex protein structures and provide the user with a good looking 3D interface in OpenGL. Based on the change of dipole moments we can construct logic circuits as described in the paper mentioned before. Furthermore I implemented several searching algorithms to find protein structures that satisfy user defined logic functions.

## Description of Branches
#### Early-searching-algorithm
Here you can find the original version of the program, as described in the [paper](https://ricsinaruto.github.io/website/docs/tdk.pdf). A very simple searching algorithm is used to find good protein structures. Basically it looks through all the possible protein structures and field value configurations until it finds one that satisfies the user defined logic function.
#### simulated-annealing
For this variant of the program I implemented a simulated annealing based searching algorithm to find field values to apply to a user defined protein structure in order to achieve a user defined logic function.
#### Annealing+old-algorithm
This variant of the program uses the previously mentioned simulated annealing method to search for field values, combined with the searching algorithm in the Early-searching-algorithm branch to automatically search for protein structures at the same time.
#### Annealing+random-structures
This variant of the program generates random protein structures and applies the simulated annealing algorithm to each structure to optimize the field values applied to the proteins. Similar to the previous branch, but instead of looking through all possible structures with increasing size, we take one random sized structure
