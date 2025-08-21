import openmc
import openmc.lib
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

##### RUN INPUTS #####

## BASIC ##

preset = 'arc2015' #'arc2015', 'marathonpaper', or 'custom'
geometry_mode = 'fullsim' #'shieldingonly' or 'reactoronly or 'fullsim'
what_to_tally = ['neutrondamage'] #'neutrondamage' and/or 'heating' are valid inputs
photons = False #include photon transport or no
weight_windows = True
batch_no = 10
particle_no = 10000
use_shield_mat = 'ti_hydride'

## ADVANCED ##

source_energy = 'mono' #'leakage' or 'mono'
damage_speed = 'fastonly' #'fastonly' means only <0.1MeV neutrons will be tracked by magnet surface neutron tallies, useful for assessing magnet damage
weight_windows_iterations = 5
ww_mesh_dim = 30
max_history_splits = 1000

##### GEOMETRY FILE #####

if geometry_mode == 'fullsim':
    suffix = "whole"
elif geometry_mode == 'shieldingonly':
    suffix = "shieldwithmagnets"
elif geometry_mode == 'reactoronly':
    suffix = "onlyreactor"
else:
    raise ValueError("Invalid run mode specified. Run mode must be of type 'fullsim', 'reactoronly', or 'shieldingonly'.")

geom_h5m_filename = f"{preset}_{suffix}.h5m"

geometry_h5m = os.path.join(results_dir, geom_h5m_filename)

print(f"Using geometry file {geometry_h5m} for run mode '{geometry_mode}'")

plasma_centre_position = 420 #centre of plasma layer

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

#channel replaced with more VV in arc 2015

inconel_comp = [("Ni", 0.525),
                ("Cr", 0.19),
                ("Nb", 0.05),
                ("Mo", 0.031),
                ("Ti", 0.009),
                ("Al", 0.005),
                ("Co", 0.005),
                ("C", 0.0004),
                ("Mn", 0.00175),
                ("Si", 0.00175),
                ("P", 0.00075),
                ("S", 0.00075),
                ("B", 0.00003),
                ("Cu", 0.0015)]

remaining_fe_frac = 1-sum(wtfrac for (element, wtfrac) in inconel_comp)
inconel = openmc.Material(name='channel_mat')
for (element, wtfrac) in inconel_comp:
    inconel.add_element(element, wtfrac, percent_type='wo')
inconel.add_element("Fe", remaining_fe_frac, percent_type='wo')
inconel.set_density('g/cm3', 8.19)

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
    raise ValueError(f"Steel alloy fractions sum to {subtotal:.5f} > 1.0!")
eurofer97_steel.add_element("Fe", remaining_percent, 'wo')
eurofer97_steel.temperature = 900.0

# tf coil material

tf_coil_mat = get_winding_material(name="tfcoil")

# Neutron shield materials
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

Pb = openmc.Material(name='shield')
Pb.add_elements_from_formula("Pb")
Pb.set_density('g/cm3', 11.34)

Ni = openmc.Material()
Ni.add_element("Ni", 1.0)
Ni.set_density('g/cm3', 8.9)

hf_hydride = openmc.Material()
hf_hydride.add_elements_from_formula("HfH2")
hf_hydride.set_density('g/cm3', 11.4)

nickel_haf = openmc.Material.mix_materials(materials=[Ni, hf_hydride],
                                           fracs=[0.6, 0.4],
                                           percent_type='vo',
                                           name='shield')

wc_haf = openmc.Material.mix_materials(materials=[WC, hf_hydride],
                                       fracs=[0.6, 0.4],
                                       percent_type='vo',
                                       name='shield')

#placeholder - effectively vacuum

placeholder = openmc.Material(name='placeholder')
placeholder.add_element("H", 1.0)
placeholder.set_density('g/cm3', 1e-12) #chosen arbitrarily but effectively 0

#full list not counting shield material
materials_list = [tungsten, vanadium_alloy, blanket_mat, placeholder, tf_coil_mat, eurofer97_steel] #basically always used

if preset == 'arc2015':
    materials_list.append(inconel)
else:
    materials_list.append(channel_mat)

if use_shield_mat == 'ti_hydride': #this is surely a bad way to do this, but its easy to read
    shield_mat = ti_hydride
elif use_shield_mat == 'zr_hydride':
    shield_mat = zr_hydride
elif use_shield_mat == 'zr_boro':
    shield_mat = zr_boro
elif use_shield_mat == 'wc':
    shield_mat = WC
elif use_shield_mat == 'pb':
    shield_mat = Pb
elif use_shield_mat == 'wc_haf':
    shield_mat = nickel_haf

materials_list.append(shield_mat)

materials = openmc.Materials(materials_list)

##### NEUTRON SOURCE #####

def make_n_ring_source(r, energy='mono', plot=False):
    n_source = openmc.IndependentSource()
    n_source.space = openmc.stats.CylindricalIndependent(
        r = openmc.stats.Discrete([r], [1.0]),
        phi = openmc.stats.Uniform(a=0, b=get_rotation_angle(deg=False)),
        z = openmc.stats.Discrete([0], [1.0])
    )
    n_source.angle = openmc.stats.Isotropic()

    if energy == 'leakage':
        n_source.energy = openmc.stats.PowerLaw(a=1e3, b=14.5e6, n=0.15) #energy^0.15 between 1 keV and 14.5 MeV
        print(f"Leakage energy neutron ring source created, radius {r}cm, for run mode '{geometry_mode}'")
    elif energy == 'mono':
        n_source.energy = openmc.stats.Discrete([14.1e6], [1.0]) #14.1MeV neutrons only
        print(f"Monoenergetic neutron ring source created, radius {r}cm, for run mode '{geometry_mode}'")
    else:
        raise ValueError("Energy spectrum for source must be either 'mono' or 'leakage'")

    #plot source as sanity check
    if plot == True:
        filename = "n_source_plotted.html"
        sourceplot = openmc_source_plotter.plot_source_position(n_source)
        sourceplot.write_html(os.path.join(results_dir, filename))
        print(f"Neutron ring source plot saved as {filename}")
    
    return n_source

def make_photon_ring_source(r):
    p_source = openmc.IndependentSource()
    p_source.space = openmc.stats.CylindricalIndependent(
        r = openmc.stats.Discrete([r], [1.0]),
        phi = openmc.stats.Uniform(a=0, b=get_rotation_angle(deg=False)),
        z = openmc.stats.Discrete([0], [1.0])
    )
    p_source.angle = openmc.stats.Isotropic()
    p_source.particle = 'photon'
    p_source.energy = openmc.stats.Discrete([2e5, 1e6, 3e6], [0.25, 0.5, 0.25])

    return p_source

n_source = make_n_ring_source(r=plasma_centre_position, energy=source_energy, plot=True) #r is inner reactor edge + half reactor thickness
if source_energy == 'leakage':
    p_leakage = 0.00079
    n_leakage = 0.00074 #both from 'reactoronly' geometry modes, tallying leakage into the shielding layer. no touchy
    total_leakage = p_leakage + n_leakage

    p_strength = p_leakage/total_leakage
    n_strength = n_leakage/total_leakage

    p_source = make_photon_ring_source(r=plasma_centre_position)
    
    p_source.strength = p_strength
    n_source.strength = n_strength

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

##### VARIANCE REDUCTION MESH #####

def make_var_mesh(dim):
    var_mesh = openmc.RegularMesh.from_domain(domain=geometry,
                                              dimension=[dim, dim, dim])
    return var_mesh
    
var_mesh = make_var_mesh(ww_mesh_dim)

if weight_windows == True:
    wwg = openmc.WeightWindowGenerator(mesh=var_mesh,
                                       energy_bounds=np.geomspace(0.02, 14.1e6, num=25),
                                       particle_type='neutron')

##### GET SURFACE IDS FOR TF COILS #####

if geometry_mode != 'reactoronly' and 'neutrondamage' in what_to_tally:
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

def bounding_cell_surface_tally(particle='neutron', name=None):
    """Returns a tally of energy-binned particles over the surface of the cell bounding the geometry"""

    if name is None:
        name = f"{particle} leakage by energy bin"
    
    cell_surfaces_dict = bounded_dag_univ.cells[10000].region.get_surfaces()
    cell_surface = [surface for surface in cell_surfaces_dict.values()][0] #should only contain one surface, so accessing it like this is fine

    surface_filter = openmc.SurfaceFilter(cell_surface)

    p_filter = openmc.ParticleFilter(particle)

    low_energy = 0.02 #eV, roughly thermal
    high_energy = 14.5e6 #eV, just above maximum energy

    e_start = np.log10(low_energy)
    e_stop = np.log10(high_energy)

    e_filter = openmc.EnergyFilter(np.logspace(start=e_start, stop=e_stop, num=50))

    surface_tally = openmc.Tally()
    surface_tally.filters = [surface_filter, p_filter, e_filter]
    surface_tally.scores = ['current']
    surface_tally.name = name

    return surface_tally

def surface_tally_from_pydagmc(surface_id, particle="neutron", damage_speed=damage_speed, name=None):
    """
    Returns an openmc.Tally object for current through a given surface
    for a specified particle type ('neutron' or 'photon').

    Params:
    -----------
    surface_id : int
        The surface ID corresponding to the surface that will be tallied over. Best found by pydagmc (open source) geometry interrogation, or Cubit (proprietary).
    particle : str
        The particle type to tally. Default is 'neutron'.
    name : str (optional)
        Name to give the tally. Defaults to '{particle} current across surface {surface_id}'.

    Returns:
    --------
    openmc.Tally
        Configured OpenMC surface tally for the requested particle.
    """

    surface_id = int(surface_id)

    surface_filter = openmc.SurfaceFilter(bins=surface_id)
    surface_filter.direction = 'both'
    p_filter = openmc.ParticleFilter(particle)

    fast_e_filter = openmc.EnergyFilter([0.02, 1e5, 14.5e6]) #bins between slow (thermal - 0.1MeV) and fast (0.1MeV-14.5MeV)

    if name is None:
        name = f"{particle} current across surface {surface_id}"

    surface_tally = openmc.Tally()
    surface_tally.filters = [surface_filter, p_filter]
    if particle=='neutron':
        if damage_speed=='fastonly':
            surface_tally.filters.append(fast_e_filter)
    surface_tally.scores = ['current']
    surface_tally.name = name

    return surface_tally

def volumetric_flux_tally_regular_mesh(mesh, id, particle="neutron", name=None):
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
    
    print(f"No. of regular mesh tally cells = {np.prod(mesh.dimension)}")
    mesh_filter = openmc.MeshFilter(mesh)
    p_filter = openmc.ParticleFilter([particle])

    if name is None:
        name = f"{particle} flux in regular mesh"

    mesh_tally = openmc.Tally(name=name)
    mesh_tally.filters = [mesh_filter, p_filter]
    mesh_tally.scores = ['flux']
    mesh_tally.id = id

    return mesh_tally

def surface_current_from_mesh(meshfile, particle="neutron", name=None):
    """
    Returns an openmc.Tally object for current over the surface of an input Unstructured Mesh"""

    raise NotImplementedError("As of 11 Aug 2025 (time of writing), openmc does not yet support unstructured meshes as the input for MeshSurfaceFilter objects")

    if name is None:
        name = f"{particle} current over mesh surface"

    dummy_sp = openmc.StatePoint(os.path.join(results_dir, meshfile))

    dummy_tally_for_mesh = dummy_sp.get_tally(scores=['flux']) #flux filter is only one present in dummy file

    dummy_mesh = dummy_tally_for_mesh.find_filter(openmc.MeshFilter).mesh

    magnet_mesh = openmc.UnstructuredMesh(filename=os.path.join(results_dir, meshfile),
                                          library='moab' #for .vtk (or .h5) files
                                          )
    surface_filter = openmc.MeshSurfaceFilter(dummy_mesh)
    p_filter = openmc.ParticleFilter(particle)

    surface_tally = openmc.Tally()
    surface_tally.filters = [surface_filter, p_filter]
    surface_tally.scores = ['current']

    return surface_tally

def volumetric_flux_from_mesh(meshfile, particle="neutron", name=None):
    """
    Returns an openmc.Tally object for flux in an input Unstructured Mesh (in .vtk format)"""


    if name is None:
        if particle=='both':
            name = f"neutron and photon flux over mesh surface"
        else:
            name = f"{particle} flux over mesh surface"

    magnet_mesh = openmc.UnstructuredMesh(filename=os.path.join(results_dir, meshfile),
                                          library='moab' #for .vtk (or .h5) files
                                          )
    mesh_filter = openmc.MeshFilter(magnet_mesh)

    if particle=='both':
        p_filter = openmc.ParticleFilter(['neutron', 'photon'])
    else:
        p_filter = openmc.ParticleFilter(particle)

    flux_tally = openmc.Tally()
    flux_tally.filters = [mesh_filter, p_filter]
    flux_tally.scores = ['flux']
    flux_tally.name = name

    return flux_tally

def heating_in_magnets(name=None):

    if name == None:
        name = f"Heating in magnet material"

    m_filter = openmc.MaterialFilter(tf_coil_mat)

    heating_tally = openmc.Tally()
    heating_tally.name = name
    heating_tally.filters = [m_filter]
    heating_tally.scores = ['heating']

    return heating_tally

tallies = openmc.Tallies()
#tally leakage energy spectra if only running reactor
if geometry_mode == 'reactoronly':
    tallies.append(bounding_cell_surface_tally(particle='neutron'))
    tallies.append(bounding_cell_surface_tally(particle='photon'))
else:
    #tally volumetric flux and surface currents
    if 'neutrondamage' in what_to_tally:
        #tallies.append(volumetric_flux_from_mesh(meshfile="magnet_mesh.vtk", particle='both'))
        for i in range(len(tf_surfaceobjs)):
            tallies.append(surface_tally_from_pydagmc(surface_id=i+1))
    #tally heating in magnet material
    if 'heating' in what_to_tally:
        tallies.append(heating_in_magnets())
if weight_windows == True:
    tallies.append(volumetric_flux_tally_regular_mesh(mesh=var_mesh,
                                                      id=69,
                                                      particle='neutron'))

for tally in tallies:
    print(f"Tally '{tally.name}' added")

##### SETTINGS #####
settings = openmc.Settings()
settings.source = [n_source]
if photons == True:
    settings.photon_transport = True
    if source_energy == 'leakage':
        settings.source.append(p_source)
if weight_windows == True:
    settings.max_history_splits = max_history_splits
    settings.weight_window_generators = wwg
settings.batches = batch_no
settings.particles = particle_no
settings.run_mode = 'fixed source'
settings.output = {'path': results_dir, 'tallies': False}  # all output files now go to results_dir

model = openmc.model.Model(geometry=geometry, settings=settings, materials=materials, tallies=tallies)
model.export_to_model_xml(path='model.xml')

if weight_windows == True:
    with openmc.lib.run_in_memory():
        regular_mesh_flux_tally = openmc.lib.tallies[69]
        wws = openmc.lib.WeightWindows.from_tally(regular_mesh_flux_tally, particle='neutron')

        for i in range(weight_windows_iterations):

            openmc.lib.run()

            wws.update_magic(regular_mesh_flux_tally)

            statepoint_name = f"statepoint_magic_{i+1}.h5"
            openmc.lib.statepoint_write(filename=os.path.join(results_dir, statepoint_name))

            openmc.lib.settings.weight_windows_on = True
else:
    model.run()
        

##### RESULTS #####

reactor_power = 1.5e9 #1500MWth
e_per_fusion = 17.6 * 1.6e-13 #17.6MeV
n_per_s = reactor_power / e_per_fusion
n_per_year = n_per_s * 60 * 60 * 24 * 365.25
#dont need to multiply by slice proportion of whole reactor since only total neutrons is affected, not neutrons per unit area

statepoint_path = os.path.join(results_dir, f"statepoint.{batch_no}.h5")
results = openmc.StatePoint(statepoint_path)

if photons == True:
    p_leakage = 0.00079
else:
    p_leakage = 0
n_leakage = 0.00074 #both from 'reactoronly' geometry modes, tallying leakage into the shielding layer. no touchy
total_leakage = p_leakage + n_leakage

p_strength = p_leakage/total_leakage
n_strength = n_leakage/total_leakage

print("Found results")

def get_area(surface_id):
    surface = tf_surfaceobjs[surface_id-1] #e.g. surface id 1 corresponds to first entry of list
    area = surface.area
    return area

def get_surface_current(surface_id, particle='neutron', per_unit_area=True, normalise=True):
    tally_name = f"{particle} current across surface {surface_id}"
    try:
        surface_tally_results = results.get_tally(name=tally_name)
        resultsdf = surface_tally_results.get_pandas_dataframe()

        if damage_speed == 'fastonly':
            fastrow = resultsdf.iloc[1] #2nd row, should be higher energy neutrons
            particle_sum = fastrow['mean']
            particle_sum_sd = fastrow['std. dev.']
        else:
            particle_sum = sum(resultsdf['mean'])
            particle_sum_sd = sum(resultsdf['std. dev.'])
        particle_sum /= n_strength #per source neutron instead of per source particle
        particle_sum_sd /= n_strength
        sd_mean_ratio = particle_sum_sd/particle_sum
        if sd_mean_ratio > 0.05:
            print(f"WARNING: Standard deviation for mean {particle} current over surface {surface_id} = {sd_mean_ratio*100}% of mean value")

        if particle == 'neutron' and normalise == False:
            particle_sum *= n_per_year

        if per_unit_area == True:
            print(f"{particle} current per cm^2 for surface {surface_id}: {particle_sum/get_area(surface_id)}")
        else:
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
        flux /= n_strength #per source neutron, not per source particle
        vtk_filename = os.path.join(results_dir, f"{particle}_flux.vtk")
        mesh.write_data_to_vtk(filename=vtk_filename, datasets={"mean": flux_1d})
        print(f"Exported {particle} flux to {vtk_filename}")
    except Exception as e:
        print(f"No {particle} mesh flux tally found or export failed: {e}")

if geometry_mode != 'reactoronly' and 'neutrondamage' in what_to_tally:
    for i in range(len(tf_surfaceobjs)):
        #print(f"Area of surface {i+1}: {get_area(i+1)}cm^2")
        get_surface_current(surface_id=i+1, per_unit_area=True)