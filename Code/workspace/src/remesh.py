import gmsh
import math
import shutil
import os

def remesh_stl(src, dst):

    gmsh.initialize()
    createGeometryAndMesh(src, dst)
    gmsh.finalize()

def createGeometryAndMesh(src, dst):
    # Clear all models and merge an STL mesh that we would like to remesh (from
    # the parent directory):

    gmsh.clear()
    gmsh.merge(src)

    # We first classify ("color") the surfaces by splitting the original surface
    # along sharp geometrical features. This will create new discrete surfaces,
    # curves and points.

    # Angle between two triangles above which an edge is considered as sharp,
    # retrieved from the ONELAB database (see below):
    # angle = gmsh.onelab.getNumber('Parameters/Angle for surface detection')[0]

    angle = 40

    # For complex geometries, patches can be too complex, too elongated or too
    # large to be parametrized; setting the following option will force the
    # creation of patches that are amenable to reparametrization:
    # forceParametrizablePatches = gmsh.onelab.getNumber(
    #     'Parameters/Create surfaces guaranteed to be parametrizable')[0]

    forceParametrizablePatches = 1

    # For open surfaces include the boundary edges in the classification
    # process:
    includeBoundary = True

    # Force curves to be split on given angle:
    curveAngle = 180
    gmsh.model.mesh.classifySurfaces(angle * math.pi / 180., includeBoundary,
                                     forceParametrizablePatches,
                                     curveAngle * math.pi / 180.)

    # Create a geometry for all the discrete curves and surfaces in the mesh, by
    # computing a parametrization for each one
    gmsh.model.mesh.createGeometry()


    # Create a volume from all the surfaces
    s = gmsh.model.getEntities(2)
    l = gmsh.model.geo.addSurfaceLoop([e[1] for e in s])
    gmsh.model.geo.addVolume([l])

    gmsh.model.geo.synchronize()

    # # We specify element sizes imposed by a size field, just because we can :-)
    # f = gmsh.model.mesh.field.add("MathEval")
    # # if gmsh.onelab.getNumber('Parameters/Apply funny mesh size field?')[0]:
    # #     gmsh.model.mesh.field.setString(f, "F", "2*Sin((x+y)/5) + 3")
    # # else:
    # #     gmsh.model.mesh.field.setString(f, "F", "4")
    #
    # gmsh.model.mesh.field.setString(f, "F", "7")
    # gmsh.model.mesh.field.setAsBackgroundMesh(f)


    #These can be varied dependant on the size of the component
    gmsh.option.setNumber("Mesh.MeshSizeMin", 3)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 7)

    gmsh.model.mesh.generate(3)
    gmsh.write(dst)

