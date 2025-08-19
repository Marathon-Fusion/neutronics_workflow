import openmc
import os
import matplotlib.pyplot as plt
from build_tokamak_with_tf_coils import get_rotation_angle
import numpy as np

batch_no = 50

print(f"Current file path: {os.path.dirname(__file__)}")
results_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'results'))
print(f"Results directory location: {results_dir}")
os.makedirs(results_dir, exist_ok=True)

statepoint_path = os.path.join(results_dir, f"statepoint.{batch_no}.h5")
results = openmc.StatePoint(statepoint_path)

print("Found results")

if results.source_present == False:
    print("No source sites present")
else:
    sources = results.source
    print(f"{len(sources)} source sites present")

##### RELATIVE SOURCE STRENGTHS #####

if results.photon_transport == True:
    p_leakage = 0.00079
n_leakage = 0.00074 #both from 'reactoronly' geometry modes, tallying leakage into the shielding layer. no touchy
total_leakage = p_leakage + n_leakage

p_strength = p_leakage/total_leakage
n_strength = n_leakage/total_leakage

##### NEUTRONS PER YEAR #####

reactor_power = 1.5e9 #1500MWth
e_per_fusion = 17.6 * 1.6e-13 #17.6MeV
n_per_s = reactor_power / e_per_fusion
n_per_year = n_per_s * 60 * 60 * 24 * 365.25
#dont need to multiply by slice proportion of whole reactor since only total neutrons is affected by the slice size, not neutrons per unit area

def get_surface_current(surface_id, particle='neutron', normalise=True):
    tally_name = f"{particle} current across surface {surface_id}"
    try:
        surface_tally_results = results.get_tally(name=tally_name)
        resultsdf = surface_tally_results.get_pandas_dataframe()
        particle_sum = sum(resultsdf['mean'])
        if normalise==False:
            particle_sum *= n_per_year
        if results.photon_transport == True:
            particle_sum /= n_strength
        print(f"{particle} current for surface {surface_id}: {particle_sum}")
    except Exception as e:
        print(f"{e}")

def get_outer_surface_leakage(particle='neutron', normalise=True):
    tally_name = f"{particle} leakage by energy bin"
    surface_leakage_tally = results.get_tally(name=tally_name)
    resultsdf = surface_leakage_tally.get_pandas_dataframe()
    lowenergies = resultsdf['energy low [eV]']
    current = resultsdf['mean']
    total_particle = sum(current)
    print(f"Total {particle} leakage rate: {total_particle}")
    plt.clf() # clears any existing plot
    plt.semilogx(lowenergies, current)
    plt.title(f"{particle} surface surrent", fontsize = 16)
    plt.xlabel(f"{particle} energy (eV)", fontsize = 14)
    plt.ylabel("Current", fontsize = 14)
    plt.tight_layout()
    plt.savefig(os.path.join(results_dir, f"{particle}_surface_current.png"))
    print(f"Saved {particle} leakage graph to {particle}_surface_current.png")

def mesh_tally_to_vtk(particle="neutron", normalise=True):
    """
    Export a mesh flux tally to VTK for the specified particle type.

    Parameters
    ----------
    particle : str, optional
        The type of particle mesh tally to export ('neutron' or 'photon').
        Default is 'neutron'.
    normalise : bool, optional
        Whether or not to normalise results per source neutron ().
        Default is True.

    Notes
    -----
    Exports a VTK file named 'neutron_flux.vtk' or 'photon_flux.vtk'.
    """

    try:
        if particle=='both':
            raise NotImplementedError("No support yet for extracting both neutron and photon flux information from the same mesh")
            tally_name = f"neutron and photon flux over mesh surface"
        else:
            tally_name = f"{particle} flux over mesh surface"
        mesh_tally_results = results.get_tally(name=tally_name) #finds the first mesh that tallies flux, needs adjusting if more than one is desired
        print("Got tally")
        mesh = mesh_tally_results.find_filter(openmc.MeshFilter).mesh
        print("Got mesh")
        flux = mesh_tally_results.get_values(scores=['flux'], value='mean')
        flux_1d = flux.squeeze() #need 1d array, not 3d
        print("Got mean flux values")
        if normalise == False:
            flux *= n_per_year
        if results.photon_transport == True:
            flux /= n_strength
        vtk_filename = os.path.join(results_dir, f"{particle}_flux.vtk")
        mesh.write_data_to_vtk(filename=vtk_filename, datasets={"mean": flux_1d})
        print(f"Exported {particle} flux to {vtk_filename}")
    except Exception as e:
        print(f"No {particle} mesh flux tally found or export failed: {e}")

def get_heating_tally(normalise=False):
    """Returns the total magnet heating in kW. If normalise is True, returns total heating in eV per source neutron"""

    tally_name = f"Heating in magnet material"

    heating_tally = results.get_tally(name=tally_name)
    resultsdf = heating_tally.get_pandas_dataframe()
    heat_per_source_particle_eV = sum(resultsdf['mean'])
    heat_per_neutron_eV = heat_per_source_particle_eV / n_strength
    heat_per_neutron_J = heat_per_neutron_eV * 1.6e-19
    if normalise == False:
        heat_W = heat_per_neutron_J * n_per_s * get_rotation_angle(deg=True)/360
        return heat_W/1000
    else:
        return heat_per_neutron_eV


# for i in range(30):
#     #print(f"Area of surface {i+1}: {get_area(i+1)}cm^2")
#     print(f"{get_surface_current(surface_id=i+1)}")

# get_outer_surface_leakage(particle='neutron')
# get_outer_surface_leakage(particle='photon')

#mesh_tally_to_vtk("neutron")
print(f"Total heating: {get_heating_tally(normalise=False)} kW")

rough_vol_estimate = 4 * 0.4 * 0.48 * 2*np.pi*3.10 #assuming magnet coils are circles with 310cm radius
vol_from_solidworks = 17786174010 * 1e-9 #output from solidworks in mm^3 for 35cm shield thickness, so converting to m^3

print(f"Total magnet volume: ~18 m^3")
print(f"Volumetric heating: ~{get_heating_tally()/18} kW/m^3")