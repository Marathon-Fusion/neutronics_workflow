import openmc
import numpy as np
import openmc_source_plotter
import os
from tape_compositions import get_winding_material
from build_tokamak_with_tf_coils import get_rotation_angle
import pydagmc

print(f"Current file path: {os.path.dirname(__file__)}")
results_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'results'))
print(f"Results directory location: {results_dir}")
os.makedirs(results_dir, exist_ok=True)

geometry_h5m = os.path.join(results_dir, "tokamak_with_tf_coils.h5m")

##### MATERIALS DEFINITION #####

#First wall#
tungsten = openmc.Material(name='tungsten')
tungsten.set_density('g/cm3', 19)
tungsten.add_element('W', 1.0)

#Structural#
vanadium_alloy = openmc.Material(name='vanadium_alloy') #VCr4Ti4 as per paper
vanadium_alloy.set_density('g/cm3', 6.05)
vanadium_alloy.add_element('V', 0.92, 'wo')
vanadium_alloy.add_element('Cr', 0.04, 'wo')
vanadium_alloy.add_element('Ti', 0.04, 'wo')

#Transmutation channel#
mercury_dens = 13.6 #liquid hg
li6_dens = 0.55 #liquid li
channel_dens = (0.85 * mercury_dens) + (0.15 * li6_dens) #85% mercury

channel_mat = openmc.Material(name='channel_mat')
channel_mat.set_density('g/cm3', channel_dens)
hg_percent = 85
# Add Hg with 90% Hg198 and 10% natural composition
hg198_percent = hg_percent * 0.9
natural_hg_percent = hg_percent * 0.1

# Add Hg198
channel_mat.add_nuclide('Hg198', hg198_percent/100.0, 'ao')

# Add natural Hg isotopes (scaled to 10% of total Hg)
# Natural abundances from ENDF/B-VII.1
natural_abundances = {
    'Hg196': 0.0015,
    'Hg198': 0.0997,
    'Hg199': 0.1687,
    'Hg200': 0.2310,
    'Hg201': 0.1318,
    'Hg202': 0.2986,
    'Hg204': 0.0687
}

# Calculate total abundance excluding Hg198
total_non_hg198_abundance = sum(abundance for isotope, abundance in natural_abundances.items() 
                                if isotope != 'Hg198')

# Scale natural abundances to the remaining 10% of Hg, excluding Hg198
for isotope, abundance in natural_abundances.items():
    if isotope != 'Hg198':  # Skip Hg198 as it's already added
        # Scale the abundance to make the sum of remaining isotopes equal to natural_hg_percent
        scaled_abundance = (abundance / total_non_hg198_abundance) * natural_hg_percent
        channel_mat.add_nuclide(isotope, scaled_abundance/100, 'ao')

#Li addition
li_percent = 15
li6_percent = li_percent * 0.9
li7_percent = li_percent * 0.1
channel_mat.add_nuclide('Li6', li6_percent/100.0, 'ao')
channel_mat.add_nuclide('Li7', li7_percent/100.0, 'ao')

#Blanket#
blanket_mat = openmc.Material(name='blanket_mat')
blanket_mat.set_density('g/cm3', 1.94)
blanket_mat.add_element('F', 4.)
blanket_mat.add_element('Be', 1.)
blanket_mat.add_nuclide('Li6', 1.8)
blanket_mat.add_nuclide('Li7', 0.2)
blanket_mat.temperature = 900.0

#Outer structure#
eurofer97_steel = openmc.Material(name='blanketouter')
eurofer97_steel.set_density('g/cm3', 7.75)
# Create composition list 
eurofer97_steel_comp_list = [("Cr", 0.0893), ("C",  0.0012), ("Mn", 0.0047),
                ("V",  0.0020), ("W",  0.0107), ("Ta", 0.0014),
                ("Ti", 0.000009), ("N",  0.00018), ("P",  0.000005),
                ("S",  0.000004), ("B",  0.000001), ("Si", 0.00006),
                ("Nb", 0.000002), ("Mo", 0.000015), ("Ni", 0.000002),
                ("Cu", 0.000003)]

subtotal = 0.0
for symbol, wt_frac in eurofer97_steel_comp_list:
    eurofer97_steel.add_element(symbol, wt_frac, 'wo')
    subtotal += wt_frac

# Make the remaining composition iron 
remaining_percent = 1.0 - subtotal
if remaining_percent < 0:
    raise ValueError(f"Alloy fractions sum to {subtotal:.5f} > 1.0!")
eurofer97_steel.add_element("Fe", remaining_percent, 'wo')
eurofer97_steel.temperature = 900.0

# tf coil material

tf_coil_mat = get_winding_material(name="tfcoil")

# Neutron shield material
ti_hydride = openmc.Material(name='shield')
ti_hydride.add_elements_from_formula("TiH2")
ti_hydride.set_density('g/cm3', 3.75) #this is density for a powder I think? not sure about packing fraction or its relevance

zr_hydride = openmc.Material(name='shield')
zr_hydride.add_elements_from_formula("ZrH2")
zr_hydride.set_density('g/cm3', 5.6)

zr_boro = openmc.Material(name='shield')
zr_boro.add_elements_from_formula("ZrB4H16")
zr_boro.set_density('g/cm3', 1.13)

WC = openmc.Material(name='shield')
WC.add_elements_from_formula("WC")
WC.set_density('g/cm3', 15.63)

#placeholder - effectively vacuum

placeholder = openmc.Material(name='placeholder')
placeholder.add_element("H", 1.0)
placeholder.set_density('g/cm3', 1e-12) #effectively 0

#full list for sim
materials_list = [tungsten, vanadium_alloy, channel_mat, blanket_mat, ti_hydride, placeholder, tf_coil_mat, eurofer97_steel]

def get_materials(materials_list):
    return openmc.Materials(materials_list)

materials = get_materials(materials_list)

##### NEUTRON SOURCE #####

def make_ring_source(r, plot=False):
    n_source = openmc.IndependentSource()
    n_source.space = openmc.stats.CylindricalIndependent(
        r = openmc.stats.Discrete([r], [1.0]),
        phi = openmc.stats.Uniform(a=0, b=get_rotation_angle(deg=False)),
        z = openmc.stats.Discrete([0], [1.0])
    )
    n_source.angle = openmc.stats.Isotropic()
    n_source.energy = openmc.stats.Discrete([14.1e6], [1.0]) #14.1MeV neutrons only

    print(f"Neutron ring source created, radius {r}cm")

    #plot source as sanity check
    if plot == True:
        filename = "n_source_plotted.html"
        sourceplot = openmc_source_plotter.plot_source_position(n_source)
        sourceplot.write_html(os.path.join(results_dir, filename))
        print(f"Neutron ring source plot saved as {filename}")
    
    return n_source

n_source = make_ring_source(r=210+440/2, plot=True) #r is inner reactor edge + half reactor thickness

##### REFLECTIVE PLANES #####

def z_rotation_matrix(angle, deg=True):
    if deg==True:
        angle *= np.pi/180 #convert to rad
    rot_mat = np.array([[np.cos(angle), -np.sin(angle), 0],
                    [np.sin(angle), np.cos(angle), 0],
                    [0, 0, 1]])
    return rot_mat

plane1_norm = np.array([[0],
                    [1],
                    [0]]) #xz plane

plane2_norm = z_rotation_matrix(angle=get_rotation_angle(deg=True), deg=True) @ plane1_norm

a1 = plane1_norm[0, 0]
b1 = plane1_norm[1, 0]
c1 = plane1_norm[2, 0]

a2 = plane2_norm[0, 0]
b2 = plane2_norm[1, 0]
c2 = plane2_norm[2, 0]

plane1 = openmc.Plane(
    a=a1,
    b=b1,
    c=c1, #plane at y=0
    d=0,
    boundary_type='reflective',
    name="plane1",
    surface_id=7001
)

plane2 = openmc.Plane(
    a=a2,
    b=b2,
    c=c2,
    d=0,
    boundary_type='reflective',
    name="plane2",
    surface_id=7002
)

# print(f"Plane 1 norm: {plane1_norm}")
# print(f"Plane 2 norm: {plane2_norm}")

##### BUILD GEOMETRY #####

dagmc_universe = openmc.DAGMCUniverse(filename=geometry_h5m, auto_geom_ids = False, universe_id = 5000)
print(f"No. of cells in DAGMC model: {dagmc_universe.n_cells}")

#returns openmc.Universe bounded by a Cell
bounded_dag_univ = dagmc_universe.bounded_universe(bounded_type='sphere', padding_distance = 1) 
#padding distance ensures reflective planes get to do their thing and the vacuum boundary of the universe doesn't eat all the neutrons
#only really necessary when using a 'box' bounded type but used anyway for robustness

border_cell_region = bounded_dag_univ.cells[10000].region & +plane1 & -plane2 #wedge contained between the two planes and the boundaries of the universe
border_cell = openmc.Cell(region=border_cell_region,
                          cell_id=6000)
border_cell.fill = bounded_dag_univ

# Creates a cell from the region and fills the cell with the dagmc geometry
geometry = openmc.Geometry([border_cell])

print("Constructed geometry")

##### GET SURFACE IDS FOR TF COILS #####

materials_model = pydagmc.Model(geometry_h5m)

tf_vols = materials_model.find_volumes_by_material('tfcoil')

#biglist is nested
tf_surfaces_biglist = []
for vol in tf_vols:
    tf_surfaces_biglist.append(vol.surfaces)

#flattens list - more efficient ways to achieve this results but whatevs
tf_surfaceobjs = []
for parentvol in tf_surfaces_biglist:
    for surface in parentvol:
        tf_surfaceobjs.append(surface)

#print(f"TF coil volumes: {tf_vols}")
#print(f"TF coil surfaces: {tf_surfaceobjs}")

print(f"No. of TF coil volumes found: {len(tf_vols)}")
print(f"No. of TF coil surfaces found: {len(tf_surfaceobjs)}")

##### TALLIES #####

def surface_tallies_from_pydagmc(surface_id, particle="neutron", name=None):
    """
    Documentation not finished
    """

    surface_id = int(surface_id)

    surface_filter = openmc.SurfaceFilter(bins=surface_id)
    surface_filter.direction = 'both'
    p_filter = openmc.ParticleFilter(particle)

    if name is None:
        name = f"{particle} current across surface {surface_id}"

    surface_tally = openmc.Tally(name=name)
    surface_tally.filters = [surface_filter, p_filter]
    surface_tally.scores = ['current']

    return surface_tally

def volumetric_flux_tally_by_grid(particle="neutron", name=None):
    """
    Returns an openmc.Tally object for flux through each element of a mesh
    for a specified particle type ('neutron' or 'photon').

    Parameters:
    -----------
    particle : str
        The type of particle to tally ('neutron' or 'photon').
    name : str (optional)
        Name to give the tally. Defaults to '{Particle} flux in mesh'.

    Returns:
    --------
    openmc.Tally
        Configured OpenMC mesh tally for the requested particle.
    """
    rgrid = np.arange(0, 945, 0.25)
    phigrid = np.linspace(0, get_rotation_angle(deg=False), num=int(get_rotation_angle(deg=True)/2)) #2 degree steps
    thetagrid = np.linspace(0, np.pi, 45) #4deg steps

    mesh = openmc.SphericalMesh(r_grid=rgrid, 
                                  phi_grid=phigrid,
                                  theta_grid=thetagrid)
    
    print(f"No. of mesh tally cells = {np.prod(mesh.dimension)}")
    mesh_filter = openmc.MeshFilter(mesh)
    p_filter = openmc.ParticleFilter([particle])

    if name is None:
        name = f"{particle.capitalize()} flux in mesh"

    mesh_tally = openmc.Tally(name=name)
    mesh_tally.filters = [mesh_filter, p_filter]
    mesh_tally.scores = ['flux']

    return mesh_tally

tallies = openmc.Tallies()
#tallies.append(volumetric_flux_tally_by_grid())
for id in range(30):
    tallies.append(surface_tallies_from_pydagmc(surface_id=id+1))

##### SETTINGS #####

settings = openmc.Settings()
#settings.photon_transport = True
settings.source = [n_source]
settings.batches = 100
settings.particles = 100000
settings.run_mode = 'fixed source'
settings.output = {'path': results_dir, 'tallies': False}  # all output files now go to results_dir

model = openmc.model.Model(geometry=geometry, settings=settings, materials=materials, tallies=tallies)
model.run()

##### RESULTS #####

statepoint_path = os.path.join(results_dir, "statepoint.100.h5")
results = openmc.StatePoint(statepoint_path)

def get_area(surface_id):
    surface = tf_surfaceobjs[surface_id-1] #e.g. surface id 1 corresponds to first entry of list
    area = surface.area
    return area

def get_surface_current(surface_id, particle='neutron', per_unit_area=True):
    tally_name = f"{particle} current across surface {surface_id}"
    try:
        surface_tally_results = results.get_tally(name=tally_name)
        resultsdf = surface_tally_results.get_pandas_dataframe()
        particle_sum = sum(resultsdf['mean'])
        if per_unit_area == True:
            print(f"{particle} current per unit area for surface {surface_id}: {particle_sum/get_area(surface_id)}")
        else:
            print(f"{particle} current for surface {surface_id}: {particle_sum}")
    except Exception as e:
        print(f"{e}")

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
        vtk_filename = os.path.join(results_dir, f"{particle}_flux.vtk")
        mesh.write_data_to_vtk(filename=vtk_filename, datasets={"mean": flux})
        print(f"Exported {particle} flux to {vtk_filename}")
    except Exception as e:
        print(f"No {particle} mesh flux tally found or export failed: {e}")

for id in range(30):
    get_surface_current(id+1, per_unit_area=False)
#mesh_tally_to_vtk("neutron")