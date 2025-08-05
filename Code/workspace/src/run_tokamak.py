import openmc
import numpy as np
import openmc_source_plotter
import os
from tape_compositions import get_winding_material
from build_tokamak_with_tf_coils import get_rotation_angle

results_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'results'))
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

tf_coil_mat = get_winding_material(name="tf coil")

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
materials_list = [tungsten, vanadium_alloy, channel_mat, blanket_mat, ti_hydride, placeholder, tf_coil_mat]

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

make_ring_source(r=210+440/2, plot=True) #r is inner reactor edge + half reactor thickness

##### BUILD GEOMETRY #####

dagmc_universe = openmc.DAGMCUniverse(geometry_h5m)
print(dagmc_universe.material_names)
bounded_dag_univ = dagmc_universe.bounded_universe()



def get_materials(materials_list):
    return openmc.Materials(materials_list)

##### CREATE TALLIES #####

# def surface_particle_current_tally(particle="neutron", name=None, layer=0):
#     """
#     Create an OpenMC Tally to measure particle current at the selected surface
#     across a range of energies.

#     Parameters
#     ----------
#     particle : str, optional
#         Particle type to tally ('neutron', 'photon', etc.). Default is 'neutron'.
#     name : str, optional
#         Optional name for the tally. If not provided, a descriptive default is used.

#     Returns
#     -------
#     openmc.Tally
#         Configured tally object for selected surface current as a function of energy.

#     Notes
#     -----
#     - The energy filter divides the spectrum into 100 logarithmic bins from
#       0.01eV up to 14.1MeV, suitable for both neutron and photon tallies.
#     - The tally reports current of the specified particle type crossing
#       the surface in both directions.
#     """
    