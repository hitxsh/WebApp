import os
from uuid import uuid4
import trimesh

import numpy as np

from .utils import analyze_particle_shape

def process_stl(input_path, output_dir, simplify_face_count=100, simplify_reduction_target=None):
    """
    Loads input_path, colors all vertices red, computes metrics, exports mesh to output_dir with unique name.
    Returns (metrics_dict, processed_filename)
    """
    mesh = trimesh.load_mesh(input_path)

    if simplify_face_count is not None:
        mesh = mesh.simplify_quadric_decimation(face_count= simplify_face_count)
    else:
        mesh = mesh.simplify_quadric_decimation(simplify_reduction_target)

    segmented_mesh, visualization_results, metrics = analyze_particle_shape(mesh, save_stl = False)

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Generate a unique filename (STL does NOT store color)
    # base = os.path.splitext(os.path.basename(input_path))[0]
    # processed_name = f"{base}_{uuid4().hex}.ply"
    processed_path = os.path.join(output_dir, "output.ply")

    segmented_mesh.export(processed_path)  # PLY will preserve colors

    return metrics


