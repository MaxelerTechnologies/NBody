# N-Body Simulation

<img src="http://appgallery.maxeler.com/v0.1/app/N-BodySimulation/icon" alt="N-Body">

## Description

This App simulates interactions between N particles under gravitational forces in space. A particle's state is described by its position(x, y, z), velocity(x, y, z), and mass. Particles can either be randomly generated or read from a CSV file. The simulation can be run for a single or multiple time- steps. After computing the acceleration applied to each particle, the velocity and position are updated. A damping factor (EPS) is introduced to prevent excessive force between two very close particles.

## Content

The repo root directory contains the following items:

- APP
- DOCS
- LICENCE.txt

### APP

Directory containing project sources.

### DOCS

Documentation of the project.
  
### LICENSE.txt

License of the project.

## Information to compile

Ensure the environment variables below are correctly set:
  * `MAXELEROSDIR`
  * `MAXCOMPILERDIR`

To compile the application, run:

    make RUNRULE="<ProfileName>"

If would like to remove the distributed maxfiles before recompiling the application run the following command before compilation:

    make RUNRULE="<ProfileName>" distclean

## Makefile targets

### build  

Compiles the application

### clean  

Removes results of compilation from build directories  

### distclean  

Removes all results of comakempilation from build directories, including all maxfiles

Brain Network on [AppGallery](http://appgallery.maxeler.com/)   

