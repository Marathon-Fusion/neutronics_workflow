import openmc
import os

batch_no = 100

print(f"Current file path: {os.path.dirname(__file__)}")
results_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'results'))
print(f"Results directory location: {results_dir}")
os.makedirs(results_dir, exist_ok=True)

##### RESULTS #####

reactor_power = 1.5e9 #1500MWth
e_per_fusion = 17.6 * 1.6e-13 #17.6MeV
n_per_s = reactor_power / e_per_fusion
n_per_year = n_per_s * 60 * 60 * 24 * 365.25
#dont need to multiply by slice proportion of whole reactor since only total neutrons is affected, not neutrons per unit area

statepoint_path = os.path.join(results_dir, f"statepoint.{batch_no}.h5")
results = openmc.StatePoint(statepoint_path)

print("Found results")

def get_surface_current(surface_id, particle='neutron'):
    tally_name = f"{particle} current across surface {surface_id}"
    try:
        surface_tally_results = results.get_tally(name=tally_name)
        resultsdf = surface_tally_results.get_pandas_dataframe()
        particle_sum = sum(resultsdf['mean'])
        print(f"{particle} current for surface {surface_id}: {particle_sum}")
    except Exception as e:
        print(f"{e}")

def mesh_tally_to_vtk(particle="neutron", normalise=True):
    """
    Export a mesh flux tally to VTK for the specified particle type.

    Parameters
    ----------
    particle : str, optional
        The type of particle mesh tally to export ('neutron' or 'photon').
        Default is 'neutron'.
    normalise : bool, optional
        Whether or not to multiply current by the number of neutrons produced per year.
        Default is False.

    Notes
    -----
    Exports a VTK file named 'neutron_flux.vtk' or 'photon_flux.vtk'.
    """

    try:
        mesh_tally_results = results.get_tally(scores=['flux'])
        print("Got tally")
        mesh = mesh_tally_results.find_filter(openmc.MeshFilter).mesh
        print("Got mesh")
        flux = mesh_tally_results.get_values(scores=['flux'], value='mean')
        flux_1d = flux.squeeze() #need 1d array, not 3d
        print("Got mean flux values")
        if normalise == False:
            flux *= n_per_year
        vtk_filename = os.path.join(results_dir, f"{particle}_flux.vtk")
        mesh.write_data_to_vtk(filename=vtk_filename, datasets={"mean": flux_1d})
        print(f"Exported {particle} flux to {vtk_filename}")
    except Exception as e:
        print(f"No {particle} mesh flux tally found or export failed: {e}")

for i in range(30):
    #print(f"Area of surface {i+1}: {get_area(i+1)}cm^2")
    print(f"{get_surface_current(surface_id=i+1)}")

mesh_tally_to_vtk("neutron")