
## What this code does

Many analyses of cosmological simulations will have to perform a step along
the lines of:
	1. read group and particle catalogs,
	2. select some or all of the groups and store some of their properties,
	3. for each group, find the particles (of some type, e.g. dark matter or gas)
           that fall within a certain distance from the group's center,
	4. process those particles' properties in some way.

Examples for this sort of algorithm include the computation of
	+ profiles (pressure, gas temperature, density, ...),
	+ cumulative group properties (integrated Compton-Y, ...)
	+ ...

The provided code lets you skip many of the implementation details, so that you can
quickly write routines to extract the data you need.

It is sufficiently flexible to accomodate a wide range of possible applications,
while providing a number of pre-implemented building blocks that accomodate for the
most common use cases.

Currently, only simulation data stored in the hdf5 format are supported.
The internal layout of the data files is assumed to be somewhat similar to the Illustris
simulations (although details can differ).


## How to use the code

Please read the
[Documentation](https://leanderthiele.github.io/group_particles/html/).
Some examples that I wrote for my own research are collected in the
[examples](examples)
directory.

In order to compile, HDF5 with C++ bindings is required.
I experienced problems when trying to compile the examples with the Intel compiler;
the GNU compiler works fine.
