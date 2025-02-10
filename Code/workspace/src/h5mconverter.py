
import cadquery as cq
import gmsh
import CAD_to_OpenMC.assembly as ab



a=ab.Assembly(['simplifiedtokamak.stp'])
a.implicit_complement = 'void'
a.run()
print(a.stl_files)

