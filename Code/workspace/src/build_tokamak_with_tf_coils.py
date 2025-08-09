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

##### RADIAL LAYERS #####

inner_reactor_edge = 210
shield_thickness = 30
inner_shield_edge = inner_reactor_edge - shield_thickness

channel_thickness = 21
blanket_inboard = 54.5 #from paper design point
blanket_outboard = 79.5 #from paper design point

radial_build=[
            (paramak.LayerType.GAP, inner_shield_edge),
            (paramak.LayerType.SOLID, 1), # placeholder
            (paramak.LayerType.GAP, 1),
            (paramak.LayerType.SOLID, shield_thickness), # neutron shield
            (paramak.LayerType.SOLID, 3), # blanketouter
            (paramak.LayerType.SOLID, blanket_inboard), # blanket inboard
            (paramak.LayerType.SOLID, 3), # structural2
            (paramak.LayerType.SOLID, channel_thickness), # channel inboard
            (paramak.LayerType.SOLID, 1), # structural1
            (paramak.LayerType.SOLID, 0.5), # first wall
            (paramak.LayerType.GAP, 6), # gap
            (paramak.LayerType.PLASMA, 240), # plasma
            (paramak.LayerType.GAP, 6), # gap
            (paramak.LayerType.SOLID, 0.5), # first wall
            (paramak.LayerType.SOLID, 1), # structural1
            (paramak.LayerType.SOLID, channel_thickness), # channel outboard
            (paramak.LayerType.SOLID, 3), #structural2
            (paramak.LayerType.SOLID, blanket_outboard), # blanket outboard
            (paramak.LayerType.SOLID, 3), # blanketouter
            (paramak.LayerType.SOLID, shield_thickness) #neutron shield
        ],

tot_reactor_thickness = sum(thickness for i, (layertype, thickness) in enumerate(radial_build[0]) if i != 0 )

##### TF COILS #####
#all distances in cm

#princeton d function tends to be unstable, be careful and always check outputs
#even after changing simple things

rotation_angle = 80

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

coil_inner_r = 120

thickness = 40 #in radial direction
gap_size = inner_reactor_edge - (coil_inner_r + thickness) #gap between inside of reactor and inner edge of magnet coil

coil_outer_r = coil_inner_r + thickness + tot_reactor_thickness + 2*gap_size - shield_thickness #designed to keep gap size the same on inboard and outboard

#again this princeton coil function is not very stable, check output is correct every time even after changing innocuous things

def main():

    print(f"Total reactor thickness (including shield): {tot_reactor_thickness}mm")
    print(f"Inner radius of TF coil: {coil_inner_r}mm")

    tf_coils = paramak.toroidal_field_coil_princeton_d(
        r1 = coil_inner_r,
        r2 = coil_outer_r,
        azimuthal_placement_angles=list(np.arange(0, rotation_angle, 20)), #20deg spacing like ARC 2015
        rotation_angle=rotation_angle,
        thickness = thickness, 
        distance=48 #half correct value
    )

    # rectangle_coil_height = 500

    # tf_coils = paramak.toroidal_field_coil_rectangle(
    #     horizontal_start_point=(coil_inner_r, rectangle_coil_height),
    #     vertical_mid_point= (coil_outer_r, 0),
    #     azimuthal_placement_angles=list(np.arange(0, rotation_angle, 20)), #20deg spacing like ARC 2015
    #     rotation_angle=rotation_angle,
    #     thickness = thickness, 
    #     distance=24 #half correct value
    # )

    print("Built TF coils...")

    ###### REACTOR #####
    #all distances in mm
    my_reactor = paramak.tokamak_from_plasma(
            radial_build=radial_build[0],
            elongation=1.6, #from paper
            triangularity=0.25, #from paper
            rotation_angle=rotation_angle, #for simplicity
            extra_cut_shapes=[tf_coils]
        )
    print("Built tokamak...")

    filename = "reactor_with_tf_coils.step"
    my_reactor.export(os.path.join(results_dir, filename))
    print(f"Tokamak model saved as {filename}")

    ##### EXPORT TO H5M #####

    A = CadToDagmc()
    #very specific order of material tags, no touchy
    A.add_cadquery_object(my_reactor, material_tags=["tfcoil", #extra_cut_shape_1
                                                    "placeholder", #layer_1
                                                    "tungsten", #layer_2
                                                    "vanadium_alloy", #layer_3
                                                    "channel_mat", #layer_4
                                                    "vanadium_alloy", #layer_5
                                                    "blanket_mat", #layer_6
                                                    "blanketouter", #layer_7
                                                    "shield", #layer_8
                                                    "placeholder"]) #plasma
    
    print("Converted CadQuery assembly to DAGMC geometry...")

    A.export_dagmc_h5m_file(filename=os.path.join(results_dir, f"tokamak_with_tf_coils.h5m"), 
                            scale_factor=1,
                            min_mesh_size=1,
                            max_mesh_size=10
                            )
    
    ##### DUPLICATE MAGNETS #####

    ids = cad_to_dagmc.get_ids_from_assembly(my_reactor)
    #print(f"IDs of assembly: {ids}")
    for_removal = [assembly for i,assembly in enumerate(ids) if i != 0] #all but the magnet assembly scheduled for removal
    print(f"For removal: {for_removal}")

    for id in for_removal:
        trimmed = my_reactor.remove(id.split("/")[-1]) #take only the assembly name after the slash
        my_reactor = copy.copy(trimmed)

    remaining_ids = cad_to_dagmc.get_ids_from_assembly(my_reactor)
    print(f"Remaining IDs: {remaining_ids}")

    B = CadToDagmc()
    B.add_cadquery_object(my_reactor, material_tags=["placeholder"]) #geometry for meshing purposes, does not need a material

    print("Created duplicate DAGMC magnet geometry...")

    B.export_unstructured_mesh_file(filename=os.path.join(results_dir, "magnet_mesh.vtk"),
                                    scale_factor=1,
                                    min_mesh_size=1,
                                    max_mesh_size=10) #same params as above, should generate the same object


if __name__ == "__main__":
    main()