import openmc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Set results directory to workspace/results
results_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'results'))
os.makedirs(results_dir, exist_ok=True)

statepoint_path = os.path.join(results_dir, "statepoint.10.h5")
results = openmc.StatePoint(statepoint_path)


def plot_surface_current(particle="neutron", normalise = True):
    """
    Plot the surface current as a function of energy for the specified particle type.
    
    Parameters
    ---------
    particle : str, optional
        The type of particle mesh tally to export ('neutron' or 'photon').
        Default is 'neutron'.
    normalise : bool, optional
        Determines whether or not to normalise results per neutron.
        Default is True.
        If False, multiplies result by the number of neutrons produced in the slice per year (calculated in magmirror).

    Notes
    -----
    Exports a .png semilog figure of the surface current. 
    """
    particle_cap = particle.capitalize()
    tally_name = f"{particle_cap} current at selected cylinder surface"
    try:
        surface_tally_results = results.get_tally(name=tally_name)
        resultsdf = surface_tally_results.get_pandas_dataframe()
        lowenergies = resultsdf["energy low [eV]"]
        if normalise == True:
            current = resultsdf["mean"]
        else:
            with open(os.path.join(results_dir, 'n_fluence_per_year.txt'), 'r') as input:
                n_fluence_per_year = input.read()
            current = resultsdf["mean"]*float(n_fluence_per_year)
        
        plt.clf() # clears any existing plot
        plt.semilogx(lowenergies, current)
        if normalise == True:
            normal_str = "normalised per source neutron"
        else:
            normal_str = "per cm^2 per year"
        plt.title(f"{particle_cap} Surface Current, {normal_str}", fontsize = 16)
        plt.xlabel(f"{particle_cap} energy (eV)", fontsize = 14)
        plt.ylabel("Current", fontsize = 14)
        plt.tight_layout()
        plt.savefig(os.path.join(results_dir, f"{particle}_surface_current.png"))
        print(f"Saved {particle} leakage graph to {particle}_surface_current.png")
    except Exception as e:
        print(f"No surface current tally for {particle}: {e}")

def sum_surface_current(particle='neutron', normalise = True):
    """
    Prints the total number of particles passing through the selected surface.
    
    Parameters
    ---------
    particle : str, optional
        The type of particle mesh tally to export ('neutron' or 'photon').
        Default is 'neutron'.
    normalise : bool, optional
        Determines whether or not to normalise results per neutron.
        Default is True.
        If False, multiplies result by the number of neutrons produced in the slice per year (calculated in magmirror).
    """
    particle_cap = particle.capitalize()
    tally_name = f"{particle_cap} current at selected cylinder surface"
    try:
        surface_tally_results = results.get_tally(name=tally_name)
        resultsdf = surface_tally_results.get_pandas_dataframe()
        particle_sum = sum(resultsdf['mean'])
        if normalise == True:
            print(f"{particle_cap} leakage rate: {particle_sum}")
        else:
            with open(os.path.join(results_dir, 'n_fluence_per_year.txt'), 'r') as input:
                n_fluence_per_year = input.read()
            particle_sum *= float(n_fluence_per_year)
            print(f"{particle_cap} fluence per year: {particle_sum} cm^-2")
    except Exception as e:
        print(f"{e}")



def mesh_tally_to_vtk(particle="neutron", normalise = True):
    """
    Export a mesh flux tally to VTK for the specified particle type.

    Parameters
    ----------
    particle : str, optional
        The type of particle mesh tally to export ('neutron' or 'photon').
        Default is 'neutron'.
    normalise : bool, optional
        Determines whether or not results are scaled by the number of neutrons emitted in a year in the relevant slice.
        Default is True, meaning results do not take this number into account.

    Notes
    -----
    Exports a VTK file named 'neutron_flux.vtk' or 'photon_flux.vtk'.
    """
    # Capitalize first letter for tally name
    particle_cap = particle.capitalize()
    tally_name = f"{particle_cap} flux in mesh"
    try:
        mesh_tally_results = results.get_tally(name=tally_name)
        mesh = mesh_tally_results.find_filter(openmc.MeshFilter).mesh
        if normalise == True:
            flux = mesh_tally_results.get_values(scores=['flux'], value='mean')
        else:
            with open(os.path.join(results_dir, 'n_fluence_per_year.txt'), 'r') as input:
                n_per_year_per_slice = input.read()
            flux = mesh_tally_results.get_values(scores=['flux'], value='mean')*float(n_per_year_per_slice)
        vtk_filename = os.path.join(results_dir, f"{particle}_flux.vtk")
        mesh.write_data_to_vtk(filename=vtk_filename, datasets={"mean": flux})
        print(f"Exported {particle} flux to {vtk_filename}")
    except Exception as e:
        print(f"No {particle} mesh flux tally found or export failed: {e}")


mesh_tally_to_vtk("neutron")
mesh_tally_to_vtk("photon")

plot_surface_current("photon", normalise=False)
plot_surface_current("neutron", normalise=False)

sum_surface_current("photon", normalise=False)
sum_surface_current("neutron", normalise=False)