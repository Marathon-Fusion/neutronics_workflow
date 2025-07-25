import openmc

# Define cylinder length
cylinder_length = 100.0

# Define the radial layers with their names and radii
radial_layers = [
    ("source + gap", 1260), #1260 thickness
    ("first wall", 1265), #5
    ("structural1", 1275), #10
    ("channel", 1485), #210, design point chosen for chrysopoeia paper
    ("structural2", 1515) #30
    ("blanket", 2015) #500, very rough estimate from design point
    ("blanketouter", 2045) #30
]

tungsten = openmc.Material(name='tungsten')
tungsten.set_density('g/cm3', 19)
tungsten.add_element('W', 1.0)

vanadium_alloy = openmc.Material(name='vanadium_alloy')
vanadium_alloy.set_density('g/cm3', 6.05)


# Axial boundaries (Z-planes)
z_min = openmc.ZPlane(z0=-cylinder_length/2, boundary_type='reflective')
z_max = openmc.ZPlane(z0=cylinder_length/2, boundary_type='reflective')

# Create cylinders for each radius
cylinders = []
for name, radius in radial_layers:
    cylinders.append(openmc.ZCylinder(r=radius))

# Create cells for each layer
cells = []
for i, (name, radius) in enumerate(radial_layers):
    if i == 0:
        # First layer (Plasma) - from center to first radius
        region = -cylinders[i] & +z_min & -z_max
    else:
        # Subsequent layers - between current cylinder and previous cylinder
        region = +cylinders[i-1] & -cylinders[i] & +z_min & -z_max
    
    # Create cell with appropriate material based on layer type
    if "gap" in name:
        # Void regions - no material (vacuum)
        cell = openmc.Cell(name=name, region=region)
    elif "first wall" in name:
        # Plasma region - no material (vacuum)
        cell = openmc.Cell(name=name, region=region)
        cell.fill()
    else:
        # Material regions - create placeholder materials
        # In a real simulation, you would define proper materials here
        cell = openmc.Cell(name=name, region=region)
        cell.fill()
    
    cells.append(cell)

# Create geometry
geometry = openmc.Geometry(cells)

# Create plots to visualize the radial structure
# XY plot (radial view)
plot_xy = geometry.plot(width=(700, 700), pixels=(1000, 1000), origin=(0, 0, 0), basis='xy')
plot_xy.figure.savefig('radial_build_xy.png')

# XZ plot (longitudinal view)
plot_xz = geometry.plot(width=(700, 700), pixels=(1000, 1000), origin=(0, 0, 0), basis='xz')
plot_xz.figure.savefig('radial_build_xz.png')

# YZ plot (longitudinal view)
plot_yz = geometry.plot(width=(700, 700), pixels=(1000, 1000), origin=(0, 0, 0), basis='yz')
plot_yz.figure.savefig('radial_build_yz.png')

# Print layer information
print("Radial build layers:")
for i, (name, radius) in enumerate(radial_layers):
    if i == 0:
        print(f"{i+1:2d}. {name:25s} - Center to {radius:6.1f} cm")
    else:
        prev_radius = radial_layers[i-1][1]
        thickness = radius - prev_radius
        print(f"{i+1:2d}. {name:25s} - {prev_radius:6.1f} to {radius:6.1f} cm (thickness: {thickness:5.2f} cm)")

print(f"\nTotal radius: {radial_layers[-1][1]:.1f} cm")
print(f"Cylinder length: {cylinder_length:.1f} cm")