# OpenMC Workflow for CAD Based Fusion Neutronics
An easy to install workflow for CAD based neutronics using OpenMC and DAGMC with a focus on Fusion. 

This workflow combines a set of open source tools to enable the creation unstructured mesh tallies in OpenMC without the need for commercial meshing software such as Coreform Cubit.

This workflow is deployed and ran in a containerised manner using Docker.

This repository contains a set of example scripts for running a full simulation from a .step file through to data processing and visualising using Paraview.

 generating .exo mesh files from .stl files said unstructured mesh tallies. A script is also included that creates .h5m files from .stp files using the tool CAD_to_OpenMC. The intermediary STL files created by the .h5m creation can be used to create .exo meshes.

 This workflow also includes a set of packages that are useful for OpenMC simulations such as Paramak, openmc_plasma_source, CAD_to_OpenMC and openmc_source_plotter.

## To Install

A pre-existing file structure is in place and installed by cloning the whole repository. The instructions below are tailored to this file structure but this is far from mandatory - handle your files as you please!

First obtain the Docker desktop app, found here:
https://www.docker.com/products/docker-desktop/

Install Paraview for data visualisation purposes, found here:
https://www.paraview.org/

Navigate to the directory you wish to install in and use as your workspace and clone the repository:

    git clone https://github.com/jakemarathon/neutronics_workflow/

Open cmd and navigate to the directory containing the docker file, from there run

    docker build -f Dockerfile -t kica/openmc-base .
This will take a while.

To set up local IDE functionality:
1. Make a virtual environment, ensure it has the packages from requirements.txt installed
2. Open cmd and navigate to the directory labled 'openmc_clone', from there run the following:
   
        git clone --recurse-submodules https://github.com/openmc-dev/openmc.git
3. Make sure your virtual environment is active!
    Navigate to the directory /openmc_clone/openmc and run the following:

        python -m pip install .

Note that the above does not allow you to run python/openmc scripts from within the IDE, this is to allow for static debugging.

To run python scripts input navigate to the 'Code' directory and input the following in cmd:

    docker run -v ${PWD}/workspace:/opt/marathon/workspace --rm kica/openmc-base src/<name_of_script>.py

## Example Usage

In the /src/ folder is a set of scripts and a .step file of a simplified example Tokamak reactor. The following instructions detail how to use these scripts to perform an OpenMC simulation tallying the damage energy & neutron flux in the different components of the Tokamak and how to process the simulation results to be visualised in Paraview.

### h5mconverter.py
This script converts a .step CAD file to a .h5m file for use as DAGMC geometry. A set of intemdiary .stl files are also created than can be used to create .exo files for use as unstructured meshes.

### mesh_converter.py
Converts .stl files to .exo files. The settings of this can be adjusted by changing the Max and Min mesh size parameters in remesh.py

### simpletokamaksim.py
Creates and runs an OpenMC neutronics simulation of the example tokamak with an example plasma source as the neutron source. Unstructured mesh tallies of all layers are taken and neutron flux and damage energy are tallied.

### tokamakadataprocessing.py
Processes the data normalising it by reactor power and by volume to output .vtk files of the neutron flux and rough estimation of DPA per full power year for visualisation in paraview.

## Notes and Quirks
Currently only python scripts are supported. Jupyter functionality will be added ASAP.

For some reason the unstructured mesh tallies made by this workflow have 'negative' volume, this is easily corrected for by multiplying the tally data by -1 as shown in the data_processing.py script. Have some suspicions as to what causes this and should hopefully be fixed soon!

Implicit volume/material tagging for DAGMC geometry can be weird. I found consistent results naming the individual body names '<material_name>_<object>' in the CAD software before exporting the entire assembly as a single .stp file
