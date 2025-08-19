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
shield_thickness = 35
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

rotation_angle = 80 #degrees, 4 coils

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

    print(f"Total reactor thickness (including shield): {tot_reactor_thickness}cm")
    print(f"Inner radius of TF coil: {coil_inner_r}cm")

    tf_coils = paramak.toroidal_field_coil_princeton_d(
        r1 = coil_inner_r,
        r2 = coil_outer_r,
        azimuthal_placement_angles=list(np.arange(0, rotation_angle, 20)), #20deg spacing like ARC 2015
        rotation_angle=rotation_angle,
        thickness = thickness, 
        distance=48
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

    stepfilename = "reactor_with_tf_coils.step"
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

        #creating list of volumes for removal
        if section=='whole':
            filename_no_ending = "tokamak_with_tf_coils"
            material_tags = ["tfcoil", #extra_cut_shape_1
                            "placeholder", #layer_1
                            "tungsten", #layer_2
                            "vanadium_alloy", #layer_3
                            "channel_mat", #layer_4
                            "vanadium_alloy", #layer_5
                            "blanket_mat", #layer_6
                            "blanketouter", #layer_7
                            "shield", #layer_8
                            "placeholder"] #plasma
        elif section=='onlyreactor':
            filename_no_ending = "tokamak_reactor_only"
            for i, assembly in enumerate(ids):
                if i == 0 or i == 8: #first volume (magnet coils) or 9th volume (shield)
                    for_removal.append(assembly)
            material_tags = ["placeholder", #layer_1
                            "tungsten", #layer_2
                            "vanadium_alloy", #layer_3
                            "channel_mat", #layer_4
                            "vanadium_alloy", #layer_5
                            "blanket_mat", #layer_6
                            "blanketouter", #layer_7
                            "placeholder"] #plasma
        elif section=='shieldwithmagnets':
            filename_no_ending = "tokamak_shield_and_magnets"
            for i, assembly in enumerate(ids):
                if i != 0 and i != 8: #not first volume (magnet coils) or 9th volume (shield)
                    for_removal.append(assembly)
            material_tags = ["tfcoil",
                             "shield"]
        elif section=='onlymagnets':
            filename_no_ending = "tokamak_magnets_only"
            for_removal = [assembly for i,assembly in enumerate(ids) if i != 0] #all but first volume
            material_tags = ["tfcoil"]
            
        #remove condemned volumes
        for id in for_removal:
            id_to_remove = id.split("/")[-1] #only volume name, not hash string in front
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

    export_tokamak_to_h5m(my_reactor, section='onlyreactor')
    export_tokamak_to_h5m(my_reactor, section='shieldwithmagnets')

if __name__ == "__main__":
    main()