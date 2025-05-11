"""
    Generates a watertight surface mesh from a corrupt surface mesh

    Parameters:
    -----------
    input_file : str 
        Path to the corrupt surface mesh 

    Returns:
    --------
    a watertight surface mesh
        All the holes that were created after the smoothing step are closed and fixed.

    """

import trimesh
import pymeshfix
import os

Path = Path

os.chdir(Path)
# Load the STL file
mesh = trimesh.load_mesh('smoothMesh.stl')

# Convert the Trimesh object to a PyMeshFix mesh
v, f = mesh.vertices, mesh.faces
meshfix = pymeshfix.MeshFix(v, f)

# Repair the mesh
meshfix.repair()

# Get the fixed mesh vertices and faces
fixed_vertices = meshfix.v
fixed_faces = meshfix.f

# Create a new Trimesh object with the fixed mesh
fixed_mesh = trimesh.Trimesh(vertices=fixed_vertices, faces=fixed_faces)

# Save the fixed mesh to a new STL file
fixed_mesh.export('repairedMesh.stl')

print("Mesh has been fixed and saved in fixed_mesh.stl'")
