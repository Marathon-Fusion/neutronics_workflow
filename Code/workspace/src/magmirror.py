import openmc
import openmc_source_plotter 
import pandas as pd
import numpy as np
import openmc_plasma_source
import os
from tape_compositions import get_winding_material

# Set results directory to workspace/results
results_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'results'))
os.makedirs(results_dir, exist_ok=True)


# Define cylinder length
cylinder_length = 100.0

mag_thickness = 40

layers_thicknesses = [
    ("source + gap", 126), 
    ("first wall", 0.5),
    ("structural1", 1),
    ("channel", 21), #design point chosen for chrysopoeia paper
    ("structural2", 3), 
    ("blanket", 50), #very rough estimate from design point of paper
    ("blanketouter", 3),
    ("magnet gap", 60), #based on almost nothing
    ("tf coil", mag_thickness) #based on ARC paper
]

outer_rad = sum([thickness for (name, thickness) in layers_thicknesses])

radial_layers = []
for i, (name, thickness) in enumerate(layers_thicknesses):
    if i == 0:
        radial_layers.append((name, thickness))
    else:
        radial_layers.append((name, thickness+radial_layers[i-1][1]))

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


# Axial boundaries for overall simulation(Z-planes)
z_min = openmc.ZPlane(z0=-cylinder_length/2, boundary_type='reflective')
z_max = openmc.ZPlane(z0=cylinder_length/2, boundary_type='reflective')

# Axial boundaries for TF coils
cable_side_length = 4 #40mm square cables from ARC paper

coil_width = 12*cable_side_length #12 cables (12x10 coil array)

TF_z_min = openmc.ZPlane(z0=-coil_width/2) #TF coil situated at z=0
TF_z_max = openmc.ZPlane(z0=coil_width/2)

# Create cylinders for each radius
cylinders = []
for i, (name, radius) in enumerate(radial_layers):
    if i != len(radial_layers)-1: #not last layer
        cylinders.append(openmc.ZCylinder(r=radius))
    else:
        cylinders.append(openmc.ZCylinder(r=radius, boundary_type='vacuum')) #last cylinder has vacuum boundary, all others are transmissive

def cell_fill(cell):
    # Create cell with appropriate material based on layer type
    if "gap" in name:
        cell.fill = None #no material, vacuum
    elif name=="first wall":
        cell.fill = tungsten
    elif "structural" in name: #applies to structural1 and structural2
        cell.fill = vanadium_alloy
    elif name=="channel":
        cell.fill = channel_mat
    elif name=="blanket":
        cell.fill = blanket_mat
    elif name=="blanketouter":
        cell.fill = eurofer97_steel
    elif name=='tf coil':
        cell.fill = get_winding_material()
    else:
        raise RuntimeError(f"Unknown material for layer '{name}'")
    return cell

# Create cells for each layer
cells = []
for i, (name, radius) in enumerate(radial_layers):
    if i == 0:
        # First layer (Plasma) - from center to first radius
        region = -cylinders[i] & +z_min & -z_max
    # elif name=='tf coil':
    #     #tf coil only extends a short width, not the whole simulation
    #     region = +cylinders[i-1] & -cylinders[i] & +TF_z_min & -TF_z_max
    #     lowervacregion = +cylinders[i-1] & -cylinders[i] & +z_min & -TF_z_min
    #     uppervacregion = +cylinders[i-1] & -cylinders[i] & +TF_z_max & z_max
    else:
        # Subsequent layers - between current cylinder and previous cylinder
        region = +cylinders[i-1] & -cylinders[i] & +z_min & -z_max
    
    cell = openmc.Cell(name=name, region=region)

    cell_fill(cell)
    
    cells.append(cell)

    # if name=="tf coil":
    #     lowervaccell = openmc.Cell(name="lowervaccell", region=lowervacregion)
    #     uppervaccell = openmc.Cell(name="uppervaccell", region=uppervacregion)
    #     cells.append(lowervaccell)
    #     cells.append(uppervaccell)

# Create geometry
geometry = openmc.Geometry(cells)

# Create plots to visualize the radial structure
# XY plot (radial view)
# plot_xy = geometry.plot(width=(2*outer_rad, 2*outer_rad), pixels=(1000, 1000), origin=(0, 0, 0), basis='xy')
# plot_xy.figure.savefig(os.path.join(results_dir, 'radial_build_xy.png'))

# # XZ plot (longitudinal view)
# plot_xz = geometry.plot(width=(500, 500), pixels=(1000, 1000), origin=(0, 0, 0), basis='xz')
# plot_xz.figure.savefig(os.path.join(results_dir, 'radial_build_xz.png'))


# # YZ plot (longitudinal view)
# plot_yz = geometry.plot(width=(500, 500), pixels=(1000, 1000), origin=(0, 0, 0), basis='yz')
# plot_yz.figure.savefig(os.path.join(results_dir, 'radial_build_yz.png'))


# Print layer information
print("Radial build layers:")
for i, (name, radius) in enumerate(radial_layers):
    if i == 0:
        print(f"{i+1:2d}. {name:25s} - Center to {radius:6.1f} cm")
    else:
        prev_radius = radial_layers[i-1][1]
        thickness = radius - prev_radius
        print(f"{i+1:2d}. {name:25s} - {prev_radius:6.1f} to {radius:6.1f} cm (thickness: {thickness:5.2f} cm)")

print(f"\nTotal radius: {radial_layers[-1][1]:.1f} cm")
print(f"Cylinder length: {cylinder_length:.1f} cm")

######################### Neutrons per second ######################

reactor_power = 1e9 #1GWth
major_r = 420 #cm, from paper design

reactor_energy = reactor_power * 365.25 * 24 * 60 * 60 #reactor (thermal) energy in J in one year
circum = 2*np.pi*major_r
slice_proportion = cylinder_length/circum #how much of full tokamak this slice represents, approx.
slice_energy = reactor_energy * slice_proportion #energy in J generated by slice in 1 year
e_per_fusion = 17.6 * 1.6e-13 #in J
n_per_year_per_slice = slice_energy/e_per_fusion

mag_r = outer_rad - mag_thickness
mag_inner_area = mag_r * cylinder_length

n_fluence_per_year = n_per_year_per_slice/mag_inner_area
print(f"Neutrons produced per year (per {cylinder_length}cm slice): {n_per_year_per_slice}")
print(f"Neutron fluence per year at magnet surface: {n_fluence_per_year} cm^-2")
n_per_year_file = os.path.join(results_dir, 'n_fluence_per_year.txt')
with open(n_per_year_file, 'w') as output:
    output.write(str(n_fluence_per_year))

######################### Neutron Source ###########################

n_source = openmc.IndependentSource()

#line source along z axis
n_source.space = openmc.stats.CartesianIndependent(openmc.stats.Discrete([0.0], [1.0]), #x=0
                                                   openmc.stats.Discrete([0.0], [1.0]), #y=0
                                                   openmc.stats.Uniform(-cylinder_length/2, cylinder_length/2)) #z up to boundaries of cylinder
n_source.angle = openmc.stats.Isotropic()
n_source.energy = openmc.stats.Discrete([14.1e6], [1.0]) #14.1MeV neutrons only
# sourceplot = plot_source_position(n_source)
# sourceplot.write_html("./n_source_plotted.html")
######################################################################

##################### Photon Source ############################
#literally the same thing as the neutron source, just replaced with photons
# photon_source = openmc.IndependentSource()
# photon_source.space = openmc.stats.CartesianIndependent(
#     openmc.stats.Discrete([0.0], [1.0]),
#     openmc.stats.Discrete([0.0], [1.0]),
#     openmc.stats.Uniform(-cylinder_length/2, cylinder_length/2)
# )
# photon_source.angle = openmc.stats.Isotropic()
# photon_source.energy = openmc.stats.Discrete([1.0e6], [1.0])  # 1MeV photons
# photon_source.particle = 'photon'
#####################################################################


########################### neutron tallies #######################################
def volumetric_flux_tally(particle="neutron", name=None):
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
    N = 360 # number of angular mesh components
    rgrid = np.arange(0, outer_rad, 0.5)
    #zgrid = np.linspace(-cylinder_length/2, cylinder_length/2, 10)
    zgrid = [-cylinder_length/2, cylinder_length/2]
    phigrid = np.linspace(0, 2*np.pi, num=N)


    mesh = openmc.CylindricalMesh(rgrid, zgrid, phigrid)
    print(f"No. of mesh cells = {np.prod(mesh.dimension)}")
    mesh_filter = openmc.MeshFilter(mesh)
    p_filter = openmc.ParticleFilter([particle])

    if name is None:
        name = f"{particle.capitalize()} flux in mesh"

    mesh_tally = openmc.Tally(name=name)
    mesh_tally.filters = [mesh_filter, p_filter]
    mesh_tally.scores = ['flux']

    return mesh_tally


def surface_particle_current_tally(particle="neutron", name=None, layer=0):
    """
    Create an OpenMC Tally to measure particle current at the outer cylindrical surface
    across a range of energies.

    Parameters
    ----------
    particle : str, optional
        Particle type to tally ('neutron', 'photon', etc.). Default is 'neutron'.
    name : str, optional
        Optional name for the tally. If not provided, a descriptive default is used.

    Returns
    -------
    openmc.Tally
        Configured tally object for outer surface current as a function of energy.

    Notes
    -----
    - The surface filter targets the outermost cylinder, corresponding to the
      vacuum boundary of the model geometry.
    - The energy filter divides the spectrum into 100 logarithmic bins from
      0.01eV up to 14.1MeV, suitable for both neutron and photon tallies.
    - The tally reports current of the specified particle type crossing
      the surface in both directions.
    """
    
    surface_filter = openmc.SurfaceFilter(bins=[cylinders[-(layer+1)]]) # Filter for the vacuum boundary at the outermost radius
    surface_filter.direction = 'both'  # could parameterize if you wish.

    # ParticleFilter for requested type (e.g., neutron or photon)
    p_filter = openmc.ParticleFilter([particle])

    # Energy filter: 100 logarithmic bins from 0.01eV up to 14.1MeV. Adjust as needed.
    e_filter = openmc.EnergyFilter(values=np.logspace(-2, 8, 100))
   
    if name is None:  # Use a descriptive default tally name if none is given
        name = f"{particle.capitalize()} current at selected cylinder surface"

    # Create and configure the tally object
    surface_tally = openmc.Tally(name=name)
    surface_tally.filters = [surface_filter, p_filter, e_filter]
    surface_tally.scores = ['current']

    return surface_tally

##########################################################################################


tallies = openmc.Tallies()
tallies.append(volumetric_flux_tally(particle="neutron")) # volumetric neutron flux tally
tallies.append(volumetric_flux_tally(particle="photon"))  # volumetric photon flux tally
tallies.append(surface_particle_current_tally(particle = "neutron", layer=1))  # magnet surface neutron flux tally
tallies.append(surface_particle_current_tally(particle = "photon", layer=1))   # magnet surface photon flux tally

settings = openmc.Settings()
settings.photon_transport = True
settings.source = [n_source]#, photon_source]
settings.batches = 10
settings.particles = 200000
settings.run_mode = 'fixed source'
settings.output = {'path': results_dir}  # all output files now go to results_dir


model = openmc.Model(geometry=geometry, settings=settings, tallies=tallies)

model.run()