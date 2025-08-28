import openmc
import numpy as np
import scipy
import os
import paramak
import cadquery
import cad_to_dagmc
import gmsh
from cad_to_dagmc import CadToDagmc
import math
import copy

# Set results directory to workspace/results
print(f"Current file path: {os.path.dirname(__file__)}")
results_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'results'))
print(f"Results directory location: {results_dir}")
os.makedirs(results_dir, exist_ok=True)

##### RUN INPUTS #####

## BASIC ##

radial_build_preset = 'custom' #'arc2015' or 'marathonpaper' or 'custom'

#plasma out
custom_layers = [(paramak.LayerType.PLASMA, 120, "placeholder"), #plasma thickness should be half actual thickness, will be doubled later
                 (paramak.LayerType.GAP, 6, None), #all other thicknesses normal
                 (paramak.LayerType.SOLID, 1, "tungsten"),
                 (paramak.LayerType.SOLID, 1, "vanadium_alloy"),
                 (paramak.LayerType.SOLID, 21, "channel_mat"),
                 (paramak.LayerType.SOLID, 3, "vanadium_alloy"),
                 (paramak.LayerType.SOLID, 55, "blanket_mat"),
                 (paramak.LayerType.SOLID, 3, "blanketouter"),
                 (paramak.LayerType.SOLID, 40, "shield")] #the 'shield' tag must stay as 'shield' for the partial geometry extraction later to work
custom_major_rad = 420

## ADVANCED ##

marathonpaper_shield_thickness = 50
coil_shield_gap = 25 #gap between inboard edge of shield and inner edge of magnet coil, 25cm is big enough to avoid overlap between tf coils and central column
tf_coil_radial_thickness = 40
tf_coil_azimuthal_thickness = 48
rotation_angle = 40 #degrees
tf_coil_placement_angle = 20 #degrees

##### MAKE RADIAL LAYERS #####

assert radial_build_preset in ('arc2015', 'marathonpaper', 'custom'), "Unrecognised radial build preset. Must be 'arc2015', 'marathonpaper', or 'custom'."

if radial_build_preset == 'marathonpaper':
    firstwall_thickness = 0.5
    channel_thickness = 21
    blanket_inboard = 54.5
    blanket_outboard = 79.5
    inner_reactor_edge = 210
    shield_thickness = marathonpaper_shield_thickness
    inner_shield_edge = inner_reactor_edge - marathonpaper_shield_thickness

if radial_build_preset == 'marathonpaper':
    radial_build=[
                (paramak.LayerType.GAP, inner_shield_edge),
                (paramak.LayerType.SOLID, 1), # placeholder
                (paramak.LayerType.GAP, 1),
                (paramak.LayerType.SOLID, shield_thickness), # neutron shield
                (paramak.LayerType.SOLID, 3), # blanketouter
                (paramak.LayerType.SOLID, blanket_inboard), # blanket inboard, from paper
                (paramak.LayerType.SOLID, 3), # structural2
                (paramak.LayerType.SOLID, channel_thickness), # channel inboard, from paper
                (paramak.LayerType.SOLID, 1), # structural1
                (paramak.LayerType.SOLID, firstwall_thickness), # first wall
                (paramak.LayerType.GAP, 6), # gap
                (paramak.LayerType.PLASMA, 240), # plasma
                (paramak.LayerType.GAP, 6), # gap
                (paramak.LayerType.SOLID, firstwall_thickness), # first wall
                (paramak.LayerType.SOLID, 1), # structural1
                (paramak.LayerType.SOLID, channel_thickness), # channel outboard, from paper
                (paramak.LayerType.SOLID, 3), #structural2
                (paramak.LayerType.SOLID, blanket_outboard), # blanket outboard, from paper
                (paramak.LayerType.SOLID, 3), # blanketouter
                (paramak.LayerType.SOLID, shield_thickness) #neutron shield
                ]
elif radial_build_preset == 'arc2015':

    channel_thickness = 20
    shield_thickness = 51
    major_rad = 420

    radial_build=[
                (paramak.LayerType.SOLID, 1), # central column
                (paramak.LayerType.GAP, 1),
                (paramak.LayerType.SOLID, shield_thickness), # neutron shield
                (paramak.LayerType.SOLID, 3), # thermal shield
                (paramak.LayerType.SOLID, 3), # blanket tank
                (paramak.LayerType.SOLID, channel_thickness), # flibe blanket
                (paramak.LayerType.SOLID, 7), # vv
                (paramak.LayerType.SOLID, 1), # first wall
                (paramak.LayerType.GAP, 3), # gap
                (paramak.LayerType.PLASMA, 214), # plasma
                (paramak.LayerType.GAP, 3), # gap
                (paramak.LayerType.SOLID, 1), # first wall
                (paramak.LayerType.SOLID, 7), # vv
                (paramak.LayerType.SOLID, channel_thickness), # flibe blanket
                (paramak.LayerType.SOLID, 3), # blanket tank
                (paramak.LayerType.SOLID, 3), # thermal shield
                (paramak.LayerType.SOLID, shield_thickness) #neutron shield
                ]
    
    halfthickness = paramak.utils.sum_up_to_plasma(radial_build)
    plasmathickness = [thickness for (type, thickness) in radial_build if type == paramak.LayerType.PLASMA][0]
    inner_shield_edge = major_rad - plasmathickness/2 - halfthickness
    radial_build.insert(0, (paramak.LayerType.GAP, inner_shield_edge))
else:
    for i, (layertype, thickness, tag) in enumerate(custom_layers):
        if i == 0:
            assert layertype == paramak.LayerType.PLASMA, "First layer of custom build must be plasma"
        if layertype == paramak.LayerType.PLASMA:
            plasmathickness = thickness*2
            radial_build = [(layertype, plasmathickness)]
        else:
            radial_build.insert(0, (layertype, thickness))
            radial_build.append((layertype, thickness))
    
    reactorhalfthickness =  paramak.utils.sum_up_to_plasma(radial_build)
    inner_shield_edge = custom_major_rad - plasmathickness/2 - reactorhalfthickness
    radial_build.insert(0, (paramak.LayerType.GAP, 1))
    radial_build.insert(0, (paramak.LayerType.SOLID, 1)) #placeholder central vertical column
    radial_build.insert(0, (paramak.LayerType.GAP, inner_shield_edge))

tot_reactor_thickness = 0
for i, (layertype, thickness) in enumerate(radial_build):
    if i != 0: #excludes inner gap
        tot_reactor_thickness += thickness

##### TF COILS #####
#all distances in cm

#princeton d function tends to be unstable, be careful and always check outputs
#even after changing simple things

def get_rotation_angle(deg = True):
    """Returns the rotation angle used.
    Parameters
    ----------
    deg : bool, optional
        Determines whether or not to use degrees. Default is 'True'.
    """
  
    if deg == True:
        return rotation_angle
    else:
        return rotation_angle*np.pi/180

#coil radii are weird
#inner radius is side closest to centre of torus on inboard side
#outer radius is ALSO side closest to centre of torus on outboard side, not side furthest away as you might expect

coil_inner_r = inner_shield_edge - tf_coil_radial_thickness - coil_shield_gap

coil_outer_r = coil_inner_r + tf_coil_radial_thickness + tot_reactor_thickness + 2*coil_shield_gap #designed to keep gap size the same on inboard and outboard

azimuthal_placement_angles = list(np.arange(0, rotation_angle, tf_coil_placement_angle))

tf_coil_number = len(azimuthal_placement_angles)

#this princeton coil function can be unstable for large dimensions (~10x current dimensions), check output is correct if 

def main():

    print(f"Total reactor thickness (including shield): {tot_reactor_thickness}cm")
    print(f"Inner radius of TF coil: {coil_inner_r}cm")

    tf_coils = paramak.toroidal_field_coil_princeton_d(
        r1 = coil_inner_r,
        r2 = coil_outer_r,
        azimuthal_placement_angles = azimuthal_placement_angles, #20deg spacing like ARC 2015
        rotation_angle = rotation_angle,
        thickness = tf_coil_radial_thickness, 
        distance = tf_coil_azimuthal_thickness
    )

    print("Built TF coils...")

    ###### REACTOR #####
    #all distances in mm
    my_reactor = paramak.tokamak_from_plasma(
            radial_build=radial_build,
            elongation=1.6, #from paper
            triangularity=0.25, #from paper
            rotation_angle=rotation_angle, #for simplicity
            extra_cut_shapes=[tf_coils]
        )
    print("Built tokamak...")

    if radial_build_preset == 'arc2015':
        stepfilename = "arc2015.step"
    elif radial_build_preset == 'marathonpaper':
        stepfilename = "marathonpaper.step"
    else:
        stepfilename = "custom_reactor.step"
    my_reactor.export(os.path.join(results_dir, stepfilename))
    print(f"Tokamak model saved as {stepfilename} for easy viewing in CAD software...")

    ##### EXPORT GEOMETRY #####
    
    def export_tokamak_to_h5m(reactor,
                              section='whole',
                              export_vtk=False,
                              scale_factor=1,
                              min_mesh_size=1,
                              max_mesh_size=10):
        
        reactor_copy = copy.copy(reactor)
        print(f"Created copy of reactor geometry for section '{section}'...")
        ids = cad_to_dagmc.get_ids_from_assembly(reactor_copy)
        
        for_removal = []

        material_tags = []

        if radial_build_preset == 'custom':
            if section == 'whole':
                for (type, thickness, tag) in custom_layers:
                    if type==paramak.LayerType.SOLID:
                        material_tags.append(tag)
                material_tags.insert(0, "placeholder") #stupid central column
                material_tags.insert(0, "tfcoil") #tf coils, always first
                material_tags.append("placeholder") #plasma, always last
            elif section == 'onlyreactor':
                raise NotImplementedError("No partial geometry exporting")
                for (type, thickness, tag) in custom_layers:
                    if type==paramak.LayerType.SOLID and tag != 'shield':
                        material_tags.append(tag)
                material_tags.append("placeholder") #plasma, always last
            elif section == 'shieldwithmagnets':
                raise NotImplementedError("No partial geometry exporting")
                for (type, thickness, tag) in custom_layers:
                    if type==paramak.LayerType.SOLID and tag == 'shield':
                        material_tags.append(tag)
                material_tags.insert(0, "tfcoil") #tf coils, always first
        elif radial_build_preset == 'arc2015':
            if section == 'whole':
                material_tags = ["placeholder", #central column
                                 "tungsten", #first wall
                                 "inconel", #vv
                                 "blanket_mat", #flibe
                                 "inconel", #blanket tank
                                 "alsiwool",
                                 "shield",
                                 "placeholder"] #plasma
                for i in range(tf_coil_number):
                    material_tags.insert(0, "tfcoil")
            else:
                raise NotImplementedError("No partial geometry exporting")
        elif radial_build_preset == 'marathonpaper': #applies to both arc2015 and marathon_paper
            if section=='whole':
                material_tags = ["tfcoil" #extra_cut_shapes_1, only one needed here since the coils tend to overlap and create only one volume
                                "placeholder", #layer_1
                                "tungsten", #layer_2
                                "vanadium_alloy", #layer_3
                                "channel_mat", #layer_4
                                "vanadium_alloy", #layer_5
                                "blanket_mat", #layer_6
                                "blanketouter", #layer_7
                                "shield", #layer_8
                                "placeholder"] #plasma
                
            #below is an example of how you could remove only certain volumes if you wanted to reimplement partial geometry extraction
            elif section=='onlyreactor':
                # for i, assembly in enumerate(ids):
                    # if i == 0 or i == 8: #first volume (magnet coils) or 9th volume (shield)
                    #     for_removal.append(assembly)
                raise NotImplementedError("No partial geometry exporting")
                material_tags = ["placeholder", #layer_1
                                "tungsten", #layer_2
                                "vanadium_alloy", #layer_3
                                "channel_mat", #layer_4
                                "vanadium_alloy", #layer_5
                                "blanket_mat", #layer_6
                                "blanketouter", #layer_7
                                "placeholder"] #plasma
            elif section=='shieldwithmagnets':
                raise NotImplementedError("No partial geometry exporting")
                for i, assembly in enumerate(ids):
                    if i != 0 and i != 8: #not first volume (magnet coils) or 9th volume (shield)
                        for_removal.append(assembly)
                material_tags = ["tfcoil",
                                "shield"]
            elif section=='onlymagnets':
                raise NotImplementedError("No partial geometry exporting")
                for_removal = [assembly for i,assembly in enumerate(ids) if i != 0] #all but first volume
                material_tags = ["tfcoil"]
            
        filename_no_ending = f"{radial_build_preset}_{section}"

        #remove condemned volumes - currently will not do anything
        for id in for_removal:
            id_to_remove = id.split("/")[-1] #only volume name, not random string in front
            trimmed_reactor = reactor_copy.remove(id_to_remove)
            print(f"Removed volume {id_to_remove} from reactor copy...")
            reactor_copy = copy.copy(trimmed_reactor)

        remaining_volumes = cad_to_dagmc.get_ids_from_assembly(reactor_copy)

        filename = f"{filename_no_ending}.h5m"
        filepath = os.path.join(results_dir, filename)

        A = CadToDagmc()
        A.add_cadquery_object(reactor_copy,
                              material_tags=material_tags
                              )
        
        print("Created DAGMC geometry from CadQuery assembly of trimmed reactor model...")

        A.export_dagmc_h5m_file(filename=filepath,
                                scale_factor=scale_factor,
                                min_mesh_size=min_mesh_size,
                                max_mesh_size=max_mesh_size)

        if export_vtk == True:
            print("Exporting .vtk mesh of trimmed reactor model...")
            vtk_filename = f"{filename_no_ending}.vtk"
            vtk_filepath = os.path.join(results_dir, vtk_filename)
            A.export_unstructured_mesh_file(filename=vtk_filepath,
                                    scale_factor=scale_factor,
                                    min_mesh_size=min_mesh_size,
                                    max_mesh_size=max_mesh_size)

    export_tokamak_to_h5m(my_reactor, section='whole', export_vtk=False)

if __name__ == "__main__":
    main()