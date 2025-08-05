import openmc
import numpy as np
import os
import paramak
import cadquery
from cad_to_dagmc import CadToDagmc
from tape_compositions import get_winding_material

# Set results directory to workspace/results
results_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'results'))
os.makedirs(results_dir, exist_ok=True)

##### RADIAL LAYERS #####

inner_reactor_edge = 2100
shield_thickness = 300
inner_shield_edge = inner_reactor_edge - shield_thickness

channel_thickness = 210
blanket_inboard = 545 #from paper design point
blanket_outboard = 795 #from paper design point

radial_build=[
            (paramak.LayerType.GAP, inner_shield_edge),
            (paramak.LayerType.SOLID, 10), # placeholder
            (paramak.LayerType.SOLID, shield_thickness), # neutron shield
            (paramak.LayerType.SOLID, 30), # blanketouter
            (paramak.LayerType.SOLID, blanket_inboard), # blanket inboard
            (paramak.LayerType.SOLID, 30), # structural2
            (paramak.LayerType.SOLID, channel_thickness), # channel inboard
            (paramak.LayerType.SOLID, 10), # structural1
            (paramak.LayerType.SOLID, 5), # first wall
            (paramak.LayerType.GAP, 60), # gap
            (paramak.LayerType.PLASMA, 2400), # plasma
            (paramak.LayerType.GAP, 60), # gap
            (paramak.LayerType.SOLID, 5), # first wall
            (paramak.LayerType.SOLID, 10), # structural1
            (paramak.LayerType.SOLID, channel_thickness), # channel outboard
            (paramak.LayerType.SOLID, 30), #structural2
            (paramak.LayerType.SOLID, blanket_outboard), # blanket outboard
            (paramak.LayerType.SOLID, 30), # blanketouter
            (paramak.LayerType.SOLID, shield_thickness) #neutron shield
        ],

tot_reactor_thickness = sum(thickness for i, (layertype, thickness) in enumerate(radial_build[0]) if i != 0 )

##### TF COILS #####
#all distances in mm

#800-1000mm inner radius is the highest value that tends to work, smaller than ideal but whatever
#EXTREMELY unstable and weird behaviour when changing this value - I recommend not ever changing this, and ALWAYS checking output in CAD or similar before using if you do have do

rotation_angle = 80.01

coil_inner_r = 800

thickness = 400 #in radial direction
gap_size = inner_reactor_edge - (coil_inner_r + thickness) #gap between inside of reactor and inner edge of magnet coil

coil_outer_r = coil_inner_r + thickness + tot_reactor_thickness + 2*gap_size - shield_thickness #designed to keep gap size the same on inboard and outboard

#again this princeton coil function is not very stable, check output is correct every time even after changing innocuous things

tf_coils = paramak.toroidal_field_coil_princeton_d(
    r1 = coil_inner_r,
    r2 = coil_outer_r,
    azimuthal_placement_angles=list(np.arange(0, rotation_angle, 20)), #20deg spacing like ARC 2015
    rotation_angle=rotation_angle,
    thickness = thickness, 
    distance=480
)

###### REACTOR #####
#all distances in mm
my_reactor = paramak.tokamak_from_plasma(
        radial_build=radial_build[0],
        elongation=1.6, #from paper
        triangularity=0.25, #from paper
        rotation_angle=rotation_angle, #for simplicity
        extra_cut_shapes=[tf_coils]
    )

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

tf_coil_mat = get_winding_material()

# Neutron shield material
ti_hydride = openmc.Material(name='TiH2')
ti_hydride.add_elements_from_formula("TiH2")
ti_hydride.set_density('g/cm3', 3.75) #this is density for a powder I think? not sure about packing fraction or its relevance

zr_hydride = openmc.Material(name='ZrH2')
zr_hydride.add_elements_from_formula("ZrH2")
zr_hydride.set_density('g/cm3', 5.6)

zr_boro = openmc.Material(name='ZrB4H16')
zr_boro.add_elements_from_formula("ZrB4H16")
zr_boro.set_density('g/cm3', 1.13)

WC = openmc.Material(name='WC')
WC.add_elements_from_formula("WC")
WC.set_density('g/cm3', 15.63)

#placeholder - effectively vacuum

placeholder = openmc.Material(name='placeholder')
placeholder.add_element("H", 1.0)
placeholder.set_density('g/cm3', 1e-12) #effectively 0

##### EXPORT TO H5M #####


A = CadToDagmc()
A.add_cadquery_object(my_reactor, material_tags=["tfcoil",
                                                "placeholder",
                                                "tungsten",
                                                "vanadium_alloy",
                                                "channel_mat",
                                                "vanadium_alloy",
                                                "blanket_mat",
                                                "blanketouter",
                                                "TiH2",
                                                "placeholder"])
A.export_dagmc_h5m_file(filename=os.path.join(results_dir, f"tokamak_with_tf_coils.h5m"), max_mesh_size=10, min_mesh_size=1, scale_factor=.1)

# ######################### Neutrons per second ######################

# reactor_power = 1e9 #1GWth
# e_per_fusion = 17.6 * 1.6e-13 #in J
# fusions_per_s = reactor_power/e_per_fusion
# n_per_s = fusions_per_s

def main():
    print(f"Total reactor thickness (including shield): {tot_reactor_thickness}mm")
    print(f"Inner radius of TF coil: {coil_inner_r}mm")
    filename = "reactor_with_tf_coils.step"
    #my_reactor.export(os.path.join(results_dir, filename))
    #print(f"Tokamak model saved as {filename}")
    print(f"{A}")

if __name__ == "main":
    main()