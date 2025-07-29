import openmc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

results = openmc.StatePoint("statepoint.10.h5")

def plot_surface_current(normalise = True):
    surface_tally_results = results.get_tally(name="Neutron current at outer cylinder surface")
    resultsdf = surface_tally_results.get_pandas_dataframe()
    #flux_tally_results = results.get_tally("Flux in last shielding cell")

    lowenergies = resultsdf["energy low [eV]"]
    if normalise == True:
        current = resultsdf["mean"]
    else:
        with open('n_per_year_per_slice.txt', 'r') as input:
            n_per_year_per_slice = input.read()
        current = resultsdf["mean"]*n_per_year_per_slice

    plt.semilogx(lowenergies, current)
    plt.xlabel("Neutron energy (eV)")
    plt.ylabel("Normalised current")
    plt.savefig("resultsoutput.png")

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
            with open('n_per_year_per_slice.txt', 'r') as input:
                n_per_year_per_slice = input.read()
            flux = mesh_tally_results.get_values(scores=['flux'], value='mean')*n_per_year_per_slice
        vtk_filename = f"{particle}_flux.vtk"
        mesh.write_data_to_vtk(filename=vtk_filename, datasets={"mean": flux})
        print(f"Exported {particle} flux to {vtk_filename}")
    except Exception as e:
        print(f"No {particle} mesh flux tally found or export failed: {e}")


mesh_tally_to_vtk("neutron")
mesh_tally_to_vtk("photon")