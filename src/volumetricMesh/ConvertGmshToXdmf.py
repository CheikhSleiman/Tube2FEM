"""
    Generates a .xdmf mesh format from a .msh mesh format

    Parameters:
    -----------
    input_file : str
        Path to the gmsh file (GmshFileName)

    Returns:
    --------
    2 Meshes 
        A surface mesh (tri.xdmf) object and a volumetric mesh (tetra.xdmf) object.

    """


import meshio
import os



def create_mesh(mesh, cell_type, prune_z=False):
    cells = mesh.get_cells_type(cell_type)

    cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
    out_mesh = meshio.Mesh(points=mesh.points, cells={cell_type: cells}, cell_data={"name_to_read":[cell_data]})
    if prune_z:
        out_mesh.prune_z_0()
    return out_mesh

# os.chdir(inputPath)
os.chdir(meshPath)
print("Current working directory:", os.getcwd()) # Confirm change

mesh3D_from_msh = meshio.read(GmshFileName)

# os.chdir(outputPath)
tetra_mesh = create_mesh(mesh3D_from_msh, "tetra")

meshio.write("Tetra.xdmf", tetra_mesh)

triangle_mesh = create_mesh(mesh3D_from_msh, "triangle", prune_z=True)
meshio.write("Tri.xdmf", triangle_mesh)

