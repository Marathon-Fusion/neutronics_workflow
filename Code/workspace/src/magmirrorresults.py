import openmc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

results = openmc.StatePoint("statepoint.10.h5")

def plot_surface_current():
    surface_tally_results = results.get_tally(name="Neutron current at outer cylinder surface")
    resultsdf = surface_tally_results.get_pandas_dataframe()
    #flux_tally_results = results.get_tally("Flux in last shielding cell")

    lowenergies = resultsdf["energy low [eV]"]
    current = resultsdf["mean"]

    plt.semilogx(lowenergies, current)
    plt.xlabel("Neutron energy (eV)")
    plt.ylabel("Normalised current")
    plt.savefig("resultsoutput.png")

def mesh_tally_to_vtk():
    mesh_tally_results = results.get_tally(name="Neutron flux in mesh")

    mesh = mesh_tally_results.find_filter(openmc.MeshFilter).mesh
    #print(mesh)
    flux = mesh_tally_results.get_values(scores=['flux'], value='mean')

    mesh.write_data_to_vtk(filename="neutron_flux.vtk", datasets={"mean": flux})

mesh_tally_to_vtk()