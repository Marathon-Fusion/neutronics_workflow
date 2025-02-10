import openmc
import math
import numpy as np
import os
from openmc_plasma_source import TokamakSource
from openmc_plasma_source.plotting import plot_tokamak_source_3D, scatter_tokamak_source
import matplotlib.pyplot as plt
import dagmc_h5m_file_inspector as di
#Materials



w = openmc.Material(name='tungsten') #first wall
w.add_nuclide('W184', 1.)
w.set_density('g/cm3', 19.250)
#w.temperature = 900.0 # temperature in Kelvin

hea1 = openmc.Material(name='structural1') #structural layers
hea1.add_element('V', 0.2,'wo')
hea1.add_element('Cr', 0.2,'wo')
hea1.add_element('W', 0.2,'wo')
hea1.add_element('Ta', 0.2,'wo')
hea1.add_element('Ti', 0.2,'wo')
hea1.set_density('g/cm3', 10.7314)
#hea.temperature = 900.0


hea2 = openmc.Material(name='structural2') #structural layers
hea2.add_element('V', 0.2,'wo')
hea2.add_element('Cr', 0.2,'wo')
hea2.add_element('W', 0.2,'wo')
hea2.add_element('Ta', 0.2,'wo')
hea2.add_element('Ti', 0.2,'wo')
hea2.set_density('g/cm3', 10.7314)
#hea.temperature = 900.0

flibe = openmc.Material(name='blanket') #breeding blanket (Li-6 enrichment 90 %)
flibe.add_element('F', 4.)
flibe.add_element('Be', 1.)
flibe.add_nuclide('Li6', 1.8)
flibe.add_nuclide('Li7', 0.2)
flibe.set_density('g/cm3', 1.94)
#flibe.temperature = 900.0

void1 = openmc.Material(name='internalspace') #hydrogen to simulate vacuum
void1.add_element('H', 1.0)
void1.set_density('g/cm3', 0.0001)

void2 = openmc.Material(name='void_comp') #hydrogen to simulate vacuum
void2.add_element('H', 1.0)
void2.set_density('g/cm3', 0.0001)

#Reactor steel
rpv_steel1 = openmc.Material(name='steel')
rpv_steel1.set_density('g/cm3', 7.9)
rpv_steel1.add_nuclide("Fe54", 0.05437098, 'wo')
rpv_steel1.add_nuclide("Fe56", 0.88500663, 'wo')
rpv_steel1.add_nuclide("Fe57", 0.0208008, 'wo')
rpv_steel1.add_nuclide("Fe58", 0.00282159, 'wo')
rpv_steel1.add_nuclide("Ni58", 0.0067198, 'wo')
rpv_steel1.add_nuclide("Ni60", 0.0026776, 'wo')
rpv_steel1.add_nuclide("Mn55", 0.01, 'wo')
rpv_steel1.add_nuclide("Cr52", 0.002092475, 'wo')
rpv_steel1.add_nuclide("C0", 0.0025, 'wo')
rpv_steel1.add_nuclide("Cu63", 0.0013696, 'wo')
#rpv_steel.temperature = 900.0


materials = openmc.Materials([w,hea1,hea2,flibe,void1,void2,rpv_steel1])
#materials.export_to_xml()


# makes use of the dagmc geometry
dag_univ = openmc.DAGMCUniverse("dagmc.h5m")

# creates an edge of universe boundary surface
vac_surf = openmc.Sphere(r=10000, surface_id=9999, boundary_type="vacuum")

# adds reflective surface for the sector model at -90 degrees

reflective_1 = openmc.Plane(
    a=math.sin(0),
    b=-math.cos(0),
    c=0.0,
    d=0.0,
    surface_id=9991,
    boundary_type="reflective",
)

# adds reflective surface for the sector model at 90 degrees
reflective_2 = openmc.Plane(
    a=math.sin(math.radians(90)),
    b=-math.cos(math.radians(90)),
    c=0.0,
    d=0.0,
    surface_id=9990,
    boundary_type="reflective",
)


# specifies the region as below the universe boundary and inside the reflective surfaces
region = -vac_surf & -reflective_1 & +reflective_2
# creates a cell from the region and fills the cell with the dagmc geometry
containing_cell = openmc.Cell(cell_id=9999, region=region, fill=dag_univ)

geometry = openmc.Geometry(root=[containing_cell])
#geometry.export_to_xml()


#Plasma Source
plasma_source = TokamakSource(
    #elongation=1.557,
    elongation=1.84,
    #ion_density_centre=1.09e20,
    ion_density_centre=1.8e20,
    ion_density_peaking_factor=1,
    ion_density_pedestal=1.05e20,
    ion_density_separatrix=1e20,
    #ion_temperature_centre=45.9,
    ion_temperature_centre=27,
    ion_temperature_peaking_factor=8.06,
    ion_temperature_pedestal=2.5,
    ion_temperature_separatrix=0.5,
    major_radius=275,
    minor_radius=80,
    pedestal_radius=0.8 * 80,
    mode="H", # 3 MODES: H, L, A. We use 'H' as suggested in [1]
    #shafranov_factor=0.44789,
    shafranov_factor=0.44789,
    triangularity=0.270,
    ion_temperature_beta=6,
    angles = (0, 0.5 * np.pi)
    )

#Settings
particle_count = 10000
batches = 5
settings = openmc.Settings()
settings.dagmc = True
settings.run_mode = 'fixed source'
settings.source = plasma_source.sources
settings.batches = batches
settings.inactive = 0
settings.particles = particle_count
#settings.export_to_xml()

tallies = openmc.Tallies()

def dpa_tally_maker(name):
    mesh = openmc.UnstructuredMesh("%s.exo" % name, library='moab')
    filter = openmc.MeshFilter(mesh)
    dpa_tally = openmc.Tally(name='%s_DPA' % name)
    dpa_tally.filters = [filter]
    dpa_tally.scores = ['damage-energy', 'flux']
    tallies.append(dpa_tally)

dpa_tally_maker('firstwall')
dpa_tally_maker('structural1')
dpa_tally_maker('internalspace')
dpa_tally_maker('structural2')
dpa_tally_maker('blankettank')
dpa_tally_maker('blanketouter')




tallies.export_to_xml()

#structural1_dpa_tally,internalspace_dpa_tally,structural2_dpa_tally,blankettank_dpa_tally,blanketouter_dpa_tally


model = openmc.model.Model(geometry, materials, settings, tallies)
model.run()
