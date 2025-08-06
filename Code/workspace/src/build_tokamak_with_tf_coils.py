import openmc
import numpy as np
import scipy
import os
import paramak
import cadquery
from cad_to_dagmc import CadToDagmc

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

coil_inner_r = 800

thickness = 400 #in radial direction
gap_size = inner_reactor_edge - (coil_inner_r + thickness) #gap between inside of reactor and inner edge of magnet coil

coil_outer_r = coil_inner_r + thickness + tot_reactor_thickness + 2*gap_size - shield_thickness #designed to keep gap size the same on inboard and outboard

#again this princeton coil function is not very stable, check output is correct every time even after changing innocuous things

def main():
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
                                                    "shield",
                                                    "placeholder"]) #very specific order, no touchy
    A.export_dagmc_h5m_file(filename=os.path.join(results_dir, f"tokamak_with_tf_coils.h5m"), max_mesh_size=15, min_mesh_size=3, scale_factor=.1)

    print(f"Total reactor thickness (including shield): {tot_reactor_thickness}mm")
    print(f"Inner radius of TF coil: {coil_inner_r}mm")
    filename = "reactor_with_tf_coils.step"
    my_reactor.export(os.path.join(results_dir, filename))
    print(f"Tokamak model saved as {filename}")

if __name__ == "__main__":
    main()