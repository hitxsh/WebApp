import trimesh
import os
from uuid import uuid4

# def process_stl(input_path, output_dir):
#     """
#     - Loads input_path
#     - Computes metrics
#     - Exports the same mesh to output_dir with a unique name
#     - Returns (metrics_dict, processed_filename)
#     """
#     mesh = trimesh.load_mesh(input_path)

#     # Compute metrics
#     metrics = {
#         "volume": mesh.volume,
#         "surface_area": mesh.area,
#         "number_of_faces": len(mesh.faces),
#         "number_of_vertices": len(mesh.vertices),
#     }

#     # Ensure output directory exists
#     os.makedirs(output_dir, exist_ok=True)

#     # Generate a unique filename
#     base = os.path.splitext(os.path.basename(input_path))[0]
#     processed_name = f"{base}_{uuid4().hex}.stl"
#     processed_path = os.path.join(output_dir, processed_name)

#     # Export mesh (STL doesn’t store colors, we’ll color it in the browser)
#     mesh.export(processed_path)

#     return metrics, processed_name

import trimesh
import os
from uuid import uuid4
import numpy as np

def process_stl(input_path, output_dir):
    """
    Loads input_path, colors all vertices red, computes metrics, exports mesh to output_dir with unique name.
    Returns (metrics_dict, processed_filename)
    """
    mesh = trimesh.load_mesh(input_path)

    # Color all vertices red (RGBA)
    red_color = np.array([255, 0, 0, 255], dtype=np.uint8)
    colors = np.tile(red_color, (len(mesh.vertices), 1))
    mesh.visual.vertex_colors = colors

    # Compute metrics
    metrics = {
        "volume": mesh.volume,
        "surface_area": mesh.area,
        "number_of_faces": len(mesh.faces),
        "number_of_vertices": len(mesh.vertices),
    }

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Generate a unique filename (STL does NOT store color)
    base = os.path.splitext(os.path.basename(input_path))[0]
    processed_name = f"{base}_{uuid4().hex}.ply"
    processed_path = os.path.join(output_dir, processed_name)

    mesh.export(processed_path)  # PLY will preserve colors

    return metrics, processed_name
