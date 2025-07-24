import openmc

# Define cylinder length
cylinder_length = 100.0

# Define the radial layers with their names and radii
radial_layers = [
    ("Plasma", 44.0),
    ("Void (halo scrape-off)", 60.0),
    ("Blanket I Tube Wall I", 60.23),
    ("Blanket I – tubular zone", 78.97),
    ("Blanket I Tube Wall II", 79.2),
    ("Void (inter-blanket gap)", 80.0),
    ("Blanket II – beam-web zone", 99.0),
    ("Void (service gap)", 101.0),
    ("Reflector", 144.0),
    ("Shield A", 154.0),
    ("Shield B", 177.0),
    ("Shield C", 185.0),
    ("Void (magnet clearance)", 192.0),
    ("Dewar (vacuum vessel)", 193.0),
    ("Kapton insulation 1", 196.0),
    ("LN₂ radiation shield", 197.0),
    ("Kapton insulation 2", 204.0),
    ("Casing", 209.0),
    ("Polyimide ground-wrap", 210.0),
    ("Winding pack", 327.3)
]

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
    if "Void" in name:
        # Void regions - no material (vacuum)
        cell = openmc.Cell(name=name, region=region)
    elif "Plasma" in name:
        # Plasma region - no material (vacuum)
        cell = openmc.Cell(name=name, region=region)
    else:
        # Material regions - create placeholder materials
        # In a real simulation, you would define proper materials here
        cell = openmc.Cell(name=name, region=region)
    
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





