import numpy as np
import openmc

#YBCO
YBCO_composition = [("Y", 0.1336), 
                    ("Ba", 0.412), 
                    ("Cu", 0.286), 
                    ("O", 0.1684)] #wt%

assert np.isclose(sum([wtfrac for (element, wtfrac) in YBCO_composition]), 1.0, atol=0.0001), "YBCO material composition fractions do not sum to 100%" #off by 0.01% is allowed

YBCO = openmc.Material(name='YBCO')
for (element, wtfrac) in YBCO_composition:
    YBCO.add_element(element, wtfrac, percent_type='wo')
YBCO.set_density('g/cm3', 6.3)

#Hastelloy

hastelloy_composition = [("Co", 0.025), 
                         ("Cr", 0.16), 
                         ("Mo", 0.16), 
                         ("Fe", 0.05), 
                         ("W", 0.04), 
                         ("Mn", 0.01), 
                         ("V", 0.0035), 
                         ("Si", 0.0008), 
                         ("C", 0.0001), 
                         ("Cu", 0.005)] #wt%, all non-nickel fractions

remaining_ni_frac = 1-sum([wtfrac for (element, wtfrac) in hastelloy_composition])

hastelloy = openmc.Material(name='Hastelloy C276')
for (element, wtfrac) in hastelloy_composition:
    hastelloy.add_element(element, wtfrac, percent_type = 'wo')
hastelloy.add_element("Ni", remaining_ni_frac, percent_type = 'wo') #fill with nickel
hastelloy.set_density('g/cm3', 8.89)

#Copper
Cu = openmc.Material(name='Copper')
Cu.add_element("Cu", 1.0)
Cu.set_density('g/cm3', 8.96)

#Silver
Ag = openmc.Material(name='Silver')
Ag.add_element("Ag", 1.0)
Ag.set_density('g/cm3', 10.49)

#overall tape
hts_tape = openmc.Material.mix_materials(
    materials=[Cu, Ag, hastelloy, YBCO],
    fracs=[0.419, 0.524, 0.017, 0.040],
    percent_type='vo'
)

#stainless 316ln
stainless_comp = [("C", 0.00017),
                  ("Mn", 0.01),
                  ("Ni", 0.12),
                  ("Cr", 0.17),
                  ("Mo", 0.025),
                  ("S", 0.00015),
                  ("Si", 0.005),
                  ("P", 0.00022),
                  ("N", 0.0013)]

remaining_fe_frac = 1-sum(wtfrac for (element, wtfrac) in stainless_comp)
stainless = openmc.Material(name='Stainless 316LN')
for (element, wtfrac) in stainless_comp:
    stainless.add_element(element, wtfrac, percent_type='wo')
stainless.add_element("Fe", remaining_fe_frac, percent_type='wo')

winding_pack = openmc.Material.mix_materials(
    materials = [stainless, Cu, hts_tape],
    fracs = [0.46, 0.46, 0.08], #roughly, from ARC paper
    percent_type='vo'
)

def get_winding_material():
    """Returns an openmc.Material object of the average composition of a TF coil winding pack according to ARC 2015"""
    return winding_pack