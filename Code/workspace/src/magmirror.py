import openmc
import openmc_source_plotter 
import pandas as pd
import numpy as np

# Define cylinder length
cylinder_length = 100.0

layers_thicknesses = [
    ("source + gap", 126), #1260 thickness
    ("first wall", 0.5),
    ("structural1", 1),
    ("channel", 21), #design point chosen for chrysopoeia paper
    ("structural2", 3), #30
    ("blanket", 50), #very rough estimate from design point of paper
    ("blanketouter", 3)
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


# Axial boundaries (Z-planes)
z_min = openmc.ZPlane(z0=-cylinder_length/2, boundary_type='reflective')
z_max = openmc.ZPlane(z0=cylinder_length/2, boundary_type='reflective')

# Create cylinders for each radius
cylinders = []
for i, (name, radius) in enumerate(radial_layers):
    if i != len(radial_layers)-1:
        cylinders.append(openmc.ZCylinder(r=radius))
    else:
        cylinders.append(openmc.ZCylinder(r=radius, boundary_type='vacuum')) #last cylinder has vacuum boundary, all others are transmissive

# Create cells for each layer
cells = []
for i, (name, radius) in enumerate(radial_layers):
    if i == 0:
        # First layer (Plasma) - from center to first radius
        region = -cylinders[i] & +z_min & -z_max
    else:
        # Subsequent layers - between current cylinder and previous cylinder
        region = +cylinders[i-1] & -cylinders[i] & +z_min & -z_max
    
    cell = openmc.Cell(name=name, region=region)
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
    else:
        raise RuntimeError(f"Unknown material for layer '{name}'")
    
    cells.append(cell)

# Create geometry
geometry = openmc.Geometry(cells)

# Create plots to visualize the radial structure
# XY plot (radial view)
plot_xy = geometry.plot(width=(5000, 5000), pixels=(1000, 1000), origin=(0, 0, 0), basis='xy')
plot_xy.figure.savefig('radial_build_xy.png')

# XZ plot (longitudinal view)
plot_xz = geometry.plot(width=(5000, 5000), pixels=(1000, 1000), origin=(0, 0, 0), basis='xz')
plot_xz.figure.savefig('radial_build_xz.png')

# YZ plot (longitudinal view)
plot_yz = geometry.plot(width=(5000, 5000), pixels=(1000, 1000), origin=(0, 0, 0), basis='yz')
plot_yz.figure.savefig('radial_build_yz.png')

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

n_source = openmc.IndependentSource()

#line source along z axis
n_source.space = openmc.stats.CartesianIndependent(openmc.stats.Discrete([0.0], [1.0]), #x=0
                                                   openmc.stats.Discrete([0.0], [1.0]), #y=0
                                                   openmc.stats.Uniform(-cylinder_length/2, cylinder_length/2)) #z up to boundaries of cylinder
n_source.angle = openmc.stats.Isotropic()
n_source.energy = openmc.stats.Discrete([14.1e6], [1.0]) #14.1MeV neutrons only
#sourceplot = plot_source_position(n_source)
#sourceplot.write_html("./n_source_plotted.html")

def volumetric_neutron_flux_tally():
    "Returns an openmc.Tally object for neutron flux through each element of a mesh"

    #set up mesh for flux visualisation
    rgrid = np.arange(0, outer_rad, 0.25)
    zgrid = np.arange(-cylinder_length/2, cylinder_length/2, 0.25)
    phigrid = np.arange(0, 2*np.pi, 1*np.pi/180) #1 degree steps

    mesh = openmc.CylindricalMesh(rgrid, zgrid, phigrid)
    mesh_filter = openmc.MeshFilter(mesh)

    n_filter = openmc.ParticleFilter(['neutron'])

    mesh_tally = openmc.Tally(name="Neutron flux in mesh")
    mesh_tally.filters = [mesh_filter, n_filter]
    mesh_tally.scores=['flux']

    return mesh_tally

def surface_neutron_current_tally():
    """Returns an openmc.Tally object for neutron current across the outer cylindrical surface with energy bins"""

    #filter around outer surface
    surface_filter = openmc.SurfaceFilter(bins=[cylinders[-1]])
    surface_filter.direction = 'both'

    #filter for neutrons
    n_filter = openmc.ParticleFilter(['neutron'])

    #filter for energies between 0.01eV (less than thermal) and 14.1MeV
    e_filter = openmc.EnergyFilter(values=np.logspace(-2, 7.15, 100))

    surface_tally = openmc.Tally(name="Neutron current at outer cylinder surface")
    surface_tally.filters = [surface_filter, n_filter, e_filter]
    surface_tally.scores = ['current']

    # flux_tally = openmc.Tally(name="Flux in last shielding cell")
    # flux_tally.filters = [openmc.CellFilter(cells[-1]), n_filter]
    # flux_tally.scores = ['flux']

    return surface_tally

tallies = openmc.Tallies()
#tallies.append(surface_neutron_current_tally())
tallies.append(volumetric_neutron_flux_tally())
#tallies.append(flux_tally)

settings = openmc.Settings()
settings.source = n_source
settings.batches = 10
settings.particles = 1000
settings.run_mode = 'fixed source'

model = openmc.Model(geometry=geometry, settings=settings, tallies=tallies)

model.run()