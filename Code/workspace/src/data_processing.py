import openmc
import matplotlib.pyplot as plt
import numpy as np
import statistics


atomic_mass_unit = 1.66054e-24 #amu in grams
fusion_power = (400e6)/4 #Example tokamak power of 400MW divided by 4 as only a quarter simulated
energy_per_fusion_reaction = 17.6e6 
eV_to_Joules = 1.60218e-19
number_of_neutrons_per_second = fusion_power / (energy_per_fusion_reaction * eV_to_Joules)
number_of_neutrons_per_year = number_of_neutrons_per_second * 60 * 60 * 24 * 365.25

tungsten_E_D = 90 # threshold displacement energy of Tungsten in eV (ASTM)
tungsten_atomic_mass_grams=3.0527348e-22
tunsten_density = 19.28
tungsten_factor = tunsten_density/(tungsten_atomic_mass_grams)

#Displacement energies for composite materials are estimates
hea_E_D = 60 #Rough Estimation for structural E D
hea_atomic_mass_grams = 103.12 * atomic_mass_unit #Rough estimate
hea_density = 10.7314
hea_factor = hea_density/(hea_atomic_mass_grams)

steel_E_D = 40 #Rough Estimation
steel_density = 7.9
steel_atomic_mass_grams = 54.94 * atomic_mass_unit
steel_factor = steel_density/(steel_atomic_mass_grams)



#Internal space considered empty so no calcualtion, Flibe is massive guess as not accurate for liquid
#Have to mulitply tally scores by -1 due to mesh volumes being negative????
sp = openmc.StatePoint("statepoint.5.h5")#Rename this to whatever the name of the statepoint file is


#Uses DPA-NRT method for calulations https://doi.org/10.1016%2F0029-5493%2875%2990035-7
#See this task in Shimwell Fusion energy workshop for better understanding: 
#https://github.com/fusion-energy/neutronics-workshop/blob/main/tasks/task_06_CSG_cell_tally_DPA/1_find_dpa.ipynb
def datagrabberdpa(name,displacement_energy,material_factor):
    dpa_tally = sp.get_tally(name="%s_DPA" % name)
    mesh = dpa_tally.find_filter(openmc.MeshFilter).mesh
    dpa_raw = (dpa_tally.get_values(scores=['damage-energy'], value='mean').reshape(mesh.dimension))*(-1)
    displacements_per_source_neutron = 0.8*dpa_raw / (2*displacement_energy)
    displacements_for_all_atoms = number_of_neutrons_per_year * displacements_per_source_neutron
    dpa_nrt = displacements_for_all_atoms / (material_factor)

    mean = (sum(dpa_nrt))/(mesh.total_volume)
    print(dpa_nrt)
    print(mean)


    

    mesh.write_data_to_vtk(
        filename="%s_DPA.vtk" % name,
        datasets={"mean": dpa_nrt}  # the first "mean" is the name of the data set label inside the vtk file
    )

datagrabberdpa('firstwall', tungsten_E_D, tungsten_factor)
datagrabberdpa('structural1', hea_E_D, hea_factor)
datagrabberdpa('structural2', hea_E_D, hea_factor)
datagrabberdpa('blanketouter', steel_E_D, steel_factor)


def datagrabberflux(name):
    flux_tally = sp.get_tally(name="%s_DPA" % name)
    mesh = flux_tally.find_filter(openmc.MeshFilter).mesh
    flux_raw = (flux_tally.get_values(scores=['flux'], value='mean').reshape(mesh.dimension))*(-1)
    flux = flux_raw*number_of_neutrons_per_second
    

   
    mesh.write_data_to_vtk(
        filename="%s_flux.vtk" % name,
        datasets={"mean": flux}  # the first "mean" is the name of the data set label inside the vtk file
    )


datagrabberflux('firstwall')
datagrabberflux('structural1')
datagrabberflux('internalspace')
datagrabberflux('structural2')
datagrabberflux('blankettank')
datagrabberflux('blanketouter')

#Use Paraview to view .vtk files!