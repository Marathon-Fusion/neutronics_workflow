import openmc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

results = openmc.StatePoint("statepoint.10.h5")


def plot_surface_current(particle="neutron"):
    """
    Plot the surface current as a function of energy for the specified particle type.
    
    Parameters
    ---------
    particle : str, optional
        The type of particle mesh tally to export ('neutron' or 'photon').
        Default is 'neutron'.

    Notes
    -----
    Exports a .png semilog figure of the surface current. 
    """
    particle_cap = particle.capitalize()
    tally_name = f"{particle_cap} current at outer cylinder surface"
    try:
        surface_tally_results = results.get_tally(name=tally_name)
        resultsdf = surface_tally_results.get_pandas_dataframe()
        lowenergies = resultsdf["energy low [eV]"]
        current = resultsdf["mean"]
        plt.clf() # clears any existing plot
        plt.semilogx(lowenergies, current)
        plt.title(f"{particle_cap} Surface Current", fontsize = 16)
        plt.xlabel(f"{particle_cap} energy (eV)", fontsize = 14)
        plt.ylabel("Normalised current", fontsize = 14)
        plt.tight_layout()
        plt.savefig(f"{particle}_surface_current.png")
    except Exception as e:
        print(f"No surface current tally for {particle}: {e}")



def mesh_tally_to_vtk(particle="neutron"):
    """
    Export a mesh flux tally to VTK for the specified particle type.

    Parameters
    ----------
    particle : str, optional
        The type of particle mesh tally to export ('neutron' or 'photon').
        Default is 'neutron'.

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
        flux = mesh_tally_results.get_values(scores=['flux'], value='mean')
        vtk_filename = f"{particle}_flux.vtk"
        mesh.write_data_to_vtk(filename=vtk_filename, datasets={"mean": flux})
        print(f"Exported {particle} flux to {vtk_filename}")
    except Exception as e:
        print(f"No {particle} mesh flux tally found or export failed: {e}")


mesh_tally_to_vtk("neutron")
mesh_tally_to_vtk("photon")

plot_surface_current("photon")
plot_surface_current("neutron")