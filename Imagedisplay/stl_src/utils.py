import numpy as np
import scipy
import trimesh
import open3d
import pyglet

import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import math
import pandas as pd
import random
import os

debugging_mode = False # Debugging mode for each step movement in centre position
#Some Utilies used in the code

def get_distance(point1, point2):
    """Calculate Euclidean distance between two points."""
    vec = point1 - point2
    return np.linalg.norm(vec)

def ensure_scalar(value):
    """Convert array to scalar if needed."""
    if hasattr(value, 'ndim') and value.ndim > 0:
        if value.size == 1:
            return value.item()
    return value


def sum_column(mat):
    return np.sum(mat, axis=0)

def normalize_vector(mat):
    epsilon = 1e-12  # small value to avoid division by zero

    if len(np.shape(mat)) == 1:
        norm = np.linalg.norm(mat)
        return mat / norm if norm > epsilon else np.zeros_like(mat)

    # mat is 2D
    squared = mat ** 2
    row_sums = np.sum(squared, axis=1)  # shape (r,)
    scaling_factors = 1.0 / np.sqrt(row_sums + epsilon)  # avoid divide by zero
    mat2 = (scaling_factors[:, np.newaxis]) * mat  # shape (r,1) * (r,c)
    return mat2


def rotate_coordinates(u_vector, v_vector, normal_vector):
    """Rotate coordinates based on normal vector."""
    rotated_u = u_vector
    rotated_v = v_vector

    epsilon = 1e-12

    normal_cross = np.cross(u_vector, v_vector) / (np.linalg.norm(np.cross(u_vector, v_vector))+epsilon)
    dot_product = np.dot(normal_vector, np.transpose(normal_cross))

    if dot_product <= -1:
        rotated_u = -rotated_u
        rotated_v = -rotated_v

    projection = normal_vector - dot_product * normal_cross
    delta_projection = (normal_cross + normal_vector) / (1 + dot_product)

    rotated_u = rotated_u - delta_projection * np.dot(projection, np.transpose(rotated_u))
    rotated_v = rotated_v - delta_projection * np.dot(projection, np.transpose(rotated_v))

    return rotated_u, rotated_v


def calculate_face_normals(faces, vertices):
    """Calculate normal vectors for each face in the mesh."""
    edge0 = vertices[faces[:, 2], :] - vertices[faces[:, 1], :]
    edge1 = vertices[faces[:, 0], :] - vertices[faces[:, 2], :]

    normals = normalize_vector(np.cross(edge0, edge1))
    return normals


def calculate_vertex_normals(faces, vertices, face_normals):
    """Calculate normal vectors for each vertex in the mesh."""
    # Edge vectors
    edge0 = np.array(vertices[faces[:, 2], :] - vertices[faces[:, 1], :])
    edge1 = np.array(vertices[faces[:, 0], :] - vertices[faces[:, 2], :])
    edge2 = np.array(vertices[faces[:, 1], :] - vertices[faces[:, 0], :])

    edge0_norm = normalize_vector(edge0)
    edge1_norm = normalize_vector(edge1)
    edge2_norm = normalize_vector(edge2)

    # Edge lengths
    len_edge0 = np.sqrt((edge0[:, 0]) ** 2 + (edge0[:, 1]) ** 2 + (edge0[:, 2]) ** 2)
    len_edge1 = np.sqrt((edge1[:, 0]) ** 2 + (edge1[:, 1]) ** 2 + (edge1[:, 2]) ** 2)
    len_edge2 = np.sqrt((edge2[:, 0]) ** 2 + (edge2[:, 1]) ** 2 + (edge2[:, 2]) ** 2)

    squared_lengths = np.array([len_edge0 ** 2, len_edge1 ** 2, len_edge2 ** 2])
    squared_lengths = np.transpose(squared_lengths)

    edge_weights = np.array([
        squared_lengths[:, 0] * (squared_lengths[:, 1] + squared_lengths[:, 2] - squared_lengths[:, 0]),
        squared_lengths[:, 1] * (squared_lengths[:, 2] + squared_lengths[:, 0] - squared_lengths[:, 1]),
        squared_lengths[:, 2] * (squared_lengths[:, 0] + squared_lengths[:, 1] - squared_lengths[:, 2])
    ])

    semi_perimeter = (len_edge0 + len_edge1 + len_edge2) / 2

    # Face area calculation using Heron's formula
    face_area = np.sqrt(semi_perimeter * (semi_perimeter - len_edge0) * (semi_perimeter - len_edge1) * (semi_perimeter - len_edge2))


    corner_area = np.zeros((np.shape(faces)[0], 3))
    vertex_area = np.zeros((np.shape(vertices)[0], 1))

    # Initialize arrays for vertex normals and coordinate systems
    vertex_normals = np.zeros((np.shape(vertices)[0], 3))
    u_coord = np.zeros((np.shape(vertices)[0], 3))
    v_coord = np.zeros((np.shape(vertices)[0], 3))

    # Calculate vertex normals and areas
    for i in range(np.shape(faces)[0]):
        # Weights for vertex normals based on face area and edge lengths
        weight_v1 = face_area[i] / ((len_edge1[i] ** 2) * (len_edge2[i] ** 2))
        weight_v2 = face_area[i] / ((len_edge0[i] ** 2) * (len_edge2[i] ** 2))
        weight_v3 = face_area[i] / ((len_edge1[i] ** 2) * (len_edge0[i] ** 2))

        # Accumulate weighted normals for each vertex
        vertex_normals[faces[i][0], :] += weight_v1 * face_normals[i, :]
        vertex_normals[faces[i][1], :] += weight_v2 * face_normals[i, :]
        vertex_normals[faces[i][2], :] += weight_v3 * face_normals[i, :]

        # Calculate corner areas based on edge weights
        if edge_weights[0][i] <= 0:
            corner_area[i][1] = -0.25 * squared_lengths[i][2] * face_area[i] / (np.dot(edge0[i, :], np.transpose(edge2[i, :])))
            corner_area[i][2] = -0.25 * squared_lengths[i][1] * face_area[i] / (np.dot(edge0[i, :], np.transpose(edge1[i, :])))
            corner_area[i][0] = face_area[i] - corner_area[i][2] - corner_area[i][1]
        elif edge_weights[1][i] <= 0:
            corner_area[i][2] = -0.25 * squared_lengths[i][0] * face_area[i] / (np.dot(edge1[i, :], np.transpose(edge0[i, :])))
            corner_area[i][0] = -0.25 * squared_lengths[i][2] * face_area[i] / (np.dot(edge1[i, :], np.transpose(edge2[i, :])))
            corner_area[i][1] = face_area[i] - corner_area[i][2] - corner_area[i][0]
        elif edge_weights[2][i] <= 0:
            corner_area[i][0] = -0.25 * squared_lengths[i][1] * face_area[i] / (np.dot(edge2[i, :], np.transpose(edge1[i, :])))
            corner_area[i][1] = -0.25 * squared_lengths[i][0] * face_area[i] / (np.dot(edge2[i, :], np.transpose(edge0[i, :])))
            corner_area[i][2] = face_area[i] - corner_area[i][1] - corner_area[i][0]
        else:
            edge_scale = 0.5 * face_area[i] / (edge_weights[0][i] + edge_weights[1][i] + edge_weights[2][i])
            corner_area[i][0] = edge_scale * (edge_weights[1][i] + edge_weights[2][i])
            corner_area[i][1] = edge_scale * (edge_weights[0][i] + edge_weights[2][i])
            corner_area[i][2] = edge_scale * (edge_weights[1][i] + edge_weights[0][i])

        # Accumulate vertex areas
        vertex_area[faces[i][0]] += corner_area[i][0]
        vertex_area[faces[i][1]] += corner_area[i][1]
        vertex_area[faces[i][2]] += corner_area[i][2]

        # Set up coordinate systems for each vertex
        u_coord[faces[i][0], :] = edge2_norm[i, :]
        u_coord[faces[i][1], :] = edge0_norm[i, :]
        u_coord[faces[i][2], :] = edge1_norm[i, :]

    # Normalize vertex normals
    vertex_normals = normalize_vector(vertex_normals)

    # Complete coordinate system setup
    epsilon = 1e-12
    for i in range(np.shape(vertices)[0]):
        u_coord[i, :] = np.cross(u_coord[i, :], vertex_normals[i, :])
        u_coord[i, :] = u_coord[i, :] / (np.linalg.norm(u_coord[i, :]) + epsilon)
        v_coord[i, :] = np.cross(vertex_normals[i, :], u_coord[i, :])

    return vertex_normals, vertex_area, corner_area, u_coord, v_coord

def project_shape_operator(u_face, v_face, normal_face, ku_orig, kuv_orig, kv_orig, u_vertex, v_vertex):
    """Project shape operator from one coordinate system to another."""
    rotated_u, rotated_v = rotate_coordinates(u_vertex, v_vertex, normal_face)

    shape_operator = np.array([[ku_orig, kuv_orig], [kuv_orig, kv_orig]])

    u1 = np.dot(rotated_u, u_face)
    v1 = np.dot(rotated_u, v_face)
    u2 = np.dot(rotated_v, u_face)
    v2 = np.dot(rotated_v, v_face)

    new_ku = np.dot(np.array([u1, v1]), np.dot(shape_operator, np.transpose(np.array([u1, v1]))))
    new_kuv = np.dot(np.array([u1, v1]), np.dot(shape_operator, np.transpose(np.array([u2, v2]))))
    new_kv = np.dot(np.array([u2, v2]), np.dot(shape_operator, np.transpose(np.array([u2, v2]))))

    return new_ku, new_kuv, new_kv

def calculate_curvature(faces, vertices, vertex_normals, face_normals, vertex_area, corner_area, u_coord, v_coord):
    """Calculate curvature metrics for the mesh."""
    # Initialize shape operators for faces and vertices
    face_shape_operators = []
    vertex_shape_operators = []

    for i in range(faces.shape[0]):
        face_shape_operators.append([[0, 0], [0, 0]])

    for i in range(vertices.shape[0]):
        vertex_shape_operators.append([[0, 0], [0, 0]])

    normal_curvature = np.zeros((1, faces.shape[0]))

    # Edge vectors
    edge0 = vertices[faces[:, 2], :] - vertices[faces[:, 1], :]
    edge1 = vertices[faces[:, 0], :] - vertices[faces[:, 2], :]
    edge2 = vertices[faces[:, 1], :] - vertices[faces[:, 0], :]

    edge0_norm = normalize_vector(edge0)

    # Weights for vertex shape operators
    vertex_face_weights = np.array(np.zeros((faces.shape[0], 3)))

    for i in range(faces.shape[0]):
        normal_face = face_normals[i, :]
        tangent = edge0_norm[i, :]
        binormal = np.cross(normal_face, tangent)
        binormal = binormal / (np.linalg.norm(binormal)) if np.linalg.norm(binormal) > 0 else np.zeros_like(binormal)

        normal0 = vertex_normals[faces[i][0], :]
        normal1 = vertex_normals[faces[i][1], :]
        normal2 = vertex_normals[faces[i][2], :]

        # System of equations for shape operator
        A = np.array([
            [np.dot(edge0[i, :], tangent), np.dot(edge0[i, :], binormal), 0],
            [0, np.dot(edge0[i, :], tangent), np.dot(edge0[i, :], binormal)],
            [np.dot(edge1[i, :], tangent), np.dot(edge1[i, :], binormal), 0],
            [0, np.dot(edge1[i, :], tangent), np.dot(edge1[i, :], binormal)],
            [np.dot(edge2[i, :], tangent), np.dot(edge2[i, :], binormal), 0],
            [0, np.dot(edge2[i, :], tangent), np.dot(edge2[i, :], binormal)]
        ])

        b = np.array([
            np.dot(normal2 - normal1, tangent),np.dot(normal2 - normal1, binormal),
            np.dot(normal0 - normal2, tangent),np.dot(normal0 - normal2, binormal),
            np.dot(normal1 - normal0, tangent),np.dot(normal1 - normal0, binormal)
        ])

        x = np.linalg.lstsq(A, b, None)

        face_shape_operators[i] = np.array([[x[0][0], x[0][1]], [x[0][1], x[0][2]]])
        face_op_dot = np.dot(np.array([1, 0]), np.dot(face_shape_operators[i], np.array([[1.], [0.]])))
        normal_curvature[0][i] = ensure_scalar(face_op_dot)
        # Calculate weights for each vertex
        vertex_face_weights[i][0] = corner_area[i][0] / ensure_scalar(vertex_area[faces[i][0]])
        vertex_face_weights[i][1] = corner_area[i][1] / ensure_scalar(vertex_area[faces[i][1]])
        vertex_face_weights[i][2] = corner_area[i][2] / ensure_scalar(vertex_area[faces[i][2]])

        # Project shape operator to each vertex's coordinate system
        for j in range(3):
            new_ku, new_kuv, new_kv = project_shape_operator(
                tangent, binormal, normal_face,
                x[0][0], x[0][1], x[0][2],
                u_coord[faces[i][j], :], v_coord[faces[i][j], :]
            )

            vertex_shape_operators[faces[i][j]] += np.dot(
                vertex_face_weights[i][j],
                np.array([[new_ku, new_kuv], [new_kuv, new_kv]])
            )

    return face_shape_operators, vertex_shape_operators, vertex_face_weights

def calculate_principal_curvatures(faces, vertices, vertex_shape_operators, u_coord, v_coord):
    """Calculate principal curvature values for each vertex."""
    principal_curvature = np.zeros((2, np.shape(vertices)[0]))

    for i in range(np.shape(vertices)[0]):
        normal_plane = np.cross(u_coord[i, :], v_coord[i, :])
        rotated_u, rotated_v = rotate_coordinates(u_coord[i, :], v_coord[i, :], normal_plane)

        ku = vertex_shape_operators[i][0][0]
        kuv = vertex_shape_operators[i][0][1]
        kv = vertex_shape_operators[i][1][1]

        cosine, sine, tangent = 1, 0, 0

        # Calculate principal directions using eigendecomposition
        if kuv != 0:
            h = 0.5 * (kv - ku) / kuv

            if h < 0:
                tangent = 1 / (h - np.sqrt(1 + h ** 2))
            else:
                tangent = 1 / (h + np.sqrt(1 + h ** 2))

            cosine = 1 / np.sqrt(1 + tangent ** 2)
            sine = tangent * cosine

        # Principal curvatures
        k1 = ku - tangent * kuv
        k2 = kv + tangent * kuv

        # Order by absolute value
        if abs(k1) <= abs(k2):
            k1, k2 = k2, k1

        principal_curvature[0][i] = k1
        principal_curvature[1][i] = k2

    return principal_curvature


def get_curvature_and_areas(faces, vertices):
    """Calculate principal curvatures and vertex areas for a mesh."""
    face_normals = calculate_face_normals(faces, vertices)
    vertex_normals, vertex_area, corner_area, u_coord, v_coord = calculate_vertex_normals(faces, vertices, face_normals)
    face_shape_ops, vertex_shape_ops, weights = calculate_curvature(
        faces, vertices, vertex_normals, face_normals, vertex_area, corner_area, u_coord, v_coord
    )
    principal_curvature = calculate_principal_curvatures(faces, vertices, vertex_shape_ops, u_coord, v_coord)

    return principal_curvature, vertex_area

def find_minimum_radius(point, vertices, min_coords, max_coords):
    """Find the minimum distance from a point to any vertex in the mesh."""
    min_radius = float('inf')

    for vertex in vertices:
        distance = get_distance(point, vertex)
        if distance < min_radius:
            min_radius = distance

    if point[0]+min_radius>max_coords[0] or point[0]-min_radius<min_coords[0] \
         or point[1]+min_radius>max_coords[1] or point[1]-min_radius<min_coords[1] \
         or point[2]+min_radius>max_coords[2] or point[2]-min_radius<min_coords[2]:
        return 0

    return min_radius


def optimize_sphere_center(center, vertices, step):
    """Optimize the center position to maximize the inscribed sphere radius."""
    min_coords = vertices.min(axis=0)
    max_coords = vertices.max(axis=0)

    # Current radius
    radius = find_minimum_radius(center, vertices, min_coords, max_coords)

    #For Debugging Purposes,
    if(debugging_mode):
        print(center, radius)

    # Step vectors
    step_x = np.array([step, 0, 0])
    step_y = np.array([0, step, 0])
    step_z = np.array([0, 0, step])

    # Try neighboring points
    radius_x_pos = 0
    radius_x_neg = 0
    radius_y_pos = 0
    radius_y_neg = 0
    radius_z_pos = 0
    radius_z_neg = 0

    if center[0] + step <= max_coords[0]:
        radius_x_pos = find_minimum_radius(center + step_x, vertices, min_coords, max_coords)

    if center[0] - step >= min_coords[0]:
        radius_x_neg = find_minimum_radius(center - step_x, vertices, min_coords, max_coords)

    if center[1] + step <= max_coords[1]:
        radius_y_pos = find_minimum_radius(center + step_y, vertices, min_coords, max_coords)

    if center[1] - step >= min_coords[1]:
        radius_y_neg = find_minimum_radius(center - step_y, vertices, min_coords, max_coords)

    if center[2] + step <= max_coords[2]:
        radius_z_pos = find_minimum_radius(center + step_z, vertices, min_coords, max_coords)

    if center[2] - step >= min_coords[2]:
        radius_z_neg = find_minimum_radius(center - step_z, vertices, min_coords, max_coords)

    # Find best direction
    best_radius = max(radius, radius_x_pos, radius_x_neg, radius_y_pos, radius_y_neg, radius_z_pos, radius_z_neg)

    if radius == best_radius:
        return center, radius

    # Recursively find better center
    if radius_x_pos == best_radius:
        new_center, new_radius = optimize_sphere_center(center + step_x, vertices, step)
        if new_radius > best_radius:
            return new_center, new_radius
        return center + step_x, best_radius

    if radius_x_neg == best_radius:
        new_center, new_radius = optimize_sphere_center(center - step_x, vertices, step)
        if new_radius > best_radius:
            return new_center, new_radius
        return center - step_x, best_radius

    if radius_y_pos == best_radius:
        new_center, new_radius = optimize_sphere_center(center + step_y, vertices, step)
        if new_radius > best_radius:
            return new_center, new_radius
        return center + step, best_radius

    if radius_y_neg == best_radius:
        new_center, new_radius = optimize_sphere_center(center - step_y, vertices, step)
        if new_radius > best_radius:
            return new_center, new_radius
        return center - step_y, best_radius

    if radius_z_pos == best_radius:
        new_center, new_radius = optimize_sphere_center(center + step_z, vertices, step)
        if new_radius > best_radius:
            return new_center, new_radius
        return center + step_z, best_radius

    if radius_z_neg == best_radius:
        new_center, new_radius = optimize_sphere_center(center - step_z, vertices, step)
        if new_radius > best_radius:
            return new_center, new_radius
        return center - step_z, best_radius

    return center, radius


def find_maximum_inscribed_sphere(vertices, step=0.1):
    """Find the largest sphere that can be inscribed in the mesh."""
    # Start at center of mesh
    center = vertices.mean(axis=0)

    # Find better center point
    best_center, best_radius = optimize_sphere_center(center, vertices, step)

    return best_center, best_radius

def classify_curvature(curvatures, vertices, inscribed_radius):
    """
    Classify vertices based on their curvature values.

    Returns dictionaries with vertex classifications based on different curvature metrics.
    """
    # Threshold for considering a point as flat
    flatness_threshold = 1 / (3 * inscribed_radius)
    simple_threshold = 0.1  # For simplified classification

    # Classification dictionaries
    max_curvature_classification = {}  # For maximum principal curvature
    mean_curvature_classification = {}  # For mean curvature
    gaussian_curvature_classification = {}  # For Gaussian curvature
    min_curvature_classification = {}  # For minimum principal curvature
    alternative_classification = {}  # Alternative classification scheme

    # Mean curvature values for later use
    mean_curvature_values = {}

    # Statistics for roundness calculation
    sum_max_curvature = 0
    max_count = 0
    sum_mean_curvature = 0
    mean_count = 0
    sum_gaussian_curvature = 0
    gaussian_count = 0
    sum_min_curvature = 0
    min_count = 0


    # Convert vertices to tuples for dictionary keys
    vertex_tuples = [(v[0], v[1], v[2]) for v in vertices]

    # Calculate curvature metrics and classify each vertex
    for i in range(len(vertices)):
        k1 = curvatures[0][i]  # Maximum principal curvature
        k2 = curvatures[1][i]  # Minimum principal curvature
        vertex = vertex_tuples[i]

        # Mean curvature (H)
        mean_curvature = (k1 + k2) / 2

        # Gaussian curvature (G)
        gaussian_curvature = k1 * k2

        # Classification based on maximum principal curvature (k1)
        max_curvature_classification[vertex] = 0
        if abs(k1) <= flatness_threshold:
            k1 = 0
        if k1 < 0:
            max_curvature_classification[vertex] = -1
        elif k1 > 0:
            max_curvature_classification[vertex] = 1
            if (1/abs(k1)) <= inscribed_radius:
                sum_max_curvature += (1/abs(k1))
                max_count += 1

        # Classification based on mean curvature (H)
        mean_curvature_classification[vertex] = 0
        if abs(mean_curvature) <= flatness_threshold:
            mean_curvature = 0
        mean_curvature_values[vertex] = mean_curvature

        if mean_curvature < 0:
            mean_curvature_classification[vertex] = -1
            if (1/abs(mean_curvature)) < inscribed_radius:
                mean_curvature_classification[vertex] = -2

        elif mean_curvature > 0:
            mean_curvature_classification[vertex] = 1
            if (1/abs(mean_curvature)) < inscribed_radius:
                mean_curvature_classification[vertex] = 2
                sum_mean_curvature += (1/abs(mean_curvature))
                mean_count += 1

        # Classification based on Gaussian curvature (G)
        gaussian_curvature_classification[vertex] = 0
        if abs(gaussian_curvature) <= flatness_threshold:
            gaussian_curvature = 0

        if gaussian_curvature < 0:
            gaussian_curvature_classification[vertex] = -1
        elif gaussian_curvature > 0:
            gaussian_curvature_classification[vertex] = 1
            if (1/abs(gaussian_curvature)) <= inscribed_radius:
                gaussian_curvature_classification[vertex] = 2
                sum_gaussian_curvature += (1/abs(gaussian_curvature))
                gaussian_count += 1

        # Classification based on minimum principal curvature (k2)
        min_curvature_classification[vertex] = 0
        if abs(k2) <= flatness_threshold:
            k2 = 0

        if k2 < 0:
            min_curvature_classification[vertex] = -1
        elif k2 > 0:
            min_curvature_classification[vertex] = 1
            if (1/abs(k2)) <= inscribed_radius:
                min_curvature_classification[vertex] = 2
                sum_min_curvature += (1/abs(k2))
                min_count += 1

        # Alternative classification
        val = 0
        if abs(k1) < simple_threshold and abs(k2) < simple_threshold:
            val = 5
        elif k1 > 0 and abs(k2) <= simple_threshold:
            val = 1
        elif k1 > 0 and k2 < 0:
            val = 2
        elif k1 > 0 and k2 > 0:
            val = 3
            if (1/k1) < inscribed_radius:
                val = 5  # Distinction between convex and corner
        elif k1 < 0 and k2 < 0:
            val = 4

        alternative_classification[vertex] = val

    # Calculate roundness metrics
    max_roundness = sum_max_curvature / (max(1, max_count) * inscribed_radius)
    min_roundness = sum_min_curvature / (max(1, min_count) * inscribed_radius)
    mean_roundness = sum_mean_curvature / (max(1, mean_count) * inscribed_radius)
    gaussian_roundness = sum_gaussian_curvature / (max(1, gaussian_count) * inscribed_radius)

    # Count distribution
    classification_counts = {0: 0, 1: 0, -1: 0}
    for k in max_curvature_classification:
        classification_counts[max_curvature_classification[k]] += 1

    return {
        'max': max_curvature_classification,
        'mean': mean_curvature_classification,
        'gaussian': gaussian_curvature_classification,
        'min': min_curvature_classification,
        'alternative': alternative_classification,
        'mean_values': mean_curvature_values,
        'roundness': {
            'max': max_roundness,
            'min': min_roundness,
            'mean': mean_roundness,
            'gaussian': gaussian_roundness
        },
        'counts': classification_counts
    }

def calculate_area_statistics(vertices, classifications, vertex_areas):
    """Calculate areas based on curvature classifications."""
    concave_area = 0
    convex_area = 0
    flat_area = 0

    # Convert vertices to tuples
    vertex_tuples = [(v[0], v[1], v[2]) for v in vertices]

    for i, vertex in enumerate(vertex_tuples):
        if classifications['mean'][vertex] > 0:
            convex_area += vertex_areas[i]
        elif classifications['mean'][vertex] < 0:
            concave_area += vertex_areas[i]
        else:
            flat_area += vertex_areas[i]

    return {
        'concave': concave_area,
        'convex': convex_area,
        'flat': flat_area
    }

def find_maximum_distance(sphere):
    """Find maximum distance between any two points in a sphere."""
    max_distance = 0
    for i in range(len(sphere.vertices)):
        for j in range(i+1, len(sphere.vertices)):
            distance = get_distance(sphere.vertices[i], sphere.vertices[j])
            if distance > max_distance:
                max_distance = distance
    return max_distance


def region_segmentation(mesh, vertices, classifications):
    """
    Segment mesh vertices by curvature classification and color-code them.
    
    Colors match the image legend:
    - Yellow: Convex regions
    - Red: Highly convex regions
    - Blue: Concave regions
    - Green: Highly concave regions
    - White: Flat regions
    
    Parameters:
    -----------
    mesh : trimesh.Trimesh
        Mesh to visualize
    vertices : array-like
        Array of vertex coordinates
    classifications : dict
        Dictionary containing classification data with 'mean' key
        
    Returns:
    --------
    mesh : trimesh.Trimesh
        Mesh with updated vertex colors
    segments : dict
        Dictionary of vertex groups by classification
    """
    segmented_mesh = mesh.copy(include_visual=True)
    # Convert vertices to tuples for dictionary keys
    vertex_tuples = [(v[0], v[1], v[2]) for v in vertices]
    
    # Initialize all vertices as white (flat regions)
    for i in range(len(vertices)):
        segmented_mesh.visual.vertex_colors[i] = [255, 255, 255, 0]  # White with full opacity
    
    # Initialize segment collections
    segments = {
        'convex': [],          # Yellow
        'highly_convex': [],   # Red
        'concave': [],         # Blue
        'highly_concave': [],  # Green
        'flat': []             # White
    }
    
    # Color vertices based on classification
    for i, vertex in enumerate(vertex_tuples):
        classification_value = classifications['mean'].get(vertex, 0)
        
        if classification_value == 0:
            # Flat regions - white (already set)
            segments['flat'].append(vertex)
            
        elif classification_value == 1:
            # Convex regions - yellow
            segmented_mesh.visual.vertex_colors[i] = [255, 255, 0, 0]
            segments['convex'].append(vertex)
            
        elif classification_value == 2:
            # Highly convex regions - red
            segmented_mesh.visual.vertex_colors[i] = [255, 0, 0, 0]
            segments['highly_convex'].append(vertex)
            
        elif classification_value == -1:
            # Concave regions - blue
            segmented_mesh.visual.vertex_colors[i] = [0, 0, 255, 0]
            segments['concave'].append(vertex)
            
        elif classification_value == -2:
            # Highly concave regions - green
            segmented_mesh.visual.vertex_colors[i] = [34, 139, 34, 0]  # Forest green
            segments['highly_concave'].append(vertex)
    
    # Print statistics
    print("=== Region Segmentation Statistics ===")
    for region_type, vertices_list in segments.items():
        print(f"{region_type.replace('_', ' ').title()}: {len(vertices_list)} vertices")
    
    return segmented_mesh, segments

    
def visualize_curvature(mesh, vertices, classifications, mean_values, inscribed_radius, inscribed_center,
                        vertex_areas, circumsphere):
    """
    Visualize mesh with colors based on curvature classification.
    Also calculates weighted metrics during visualization.
    """
    # Reset mesh colors
    for i in range(len(vertices)):
        mesh.visual.vertex_colors[i] = [255, 255, 255, 0]

    # Convert vertices to tuples
    vertex_tuples = [(v[0], v[1], v[2]) for v in vertices]

    # Calculate maximum distance in circumsphere
    max_distance = find_maximum_distance(circumsphere)

    # Initialize weighted metrics
    weighted_concave = 0
    weighted_convex = 0
    iterative_term1 = 0
    iterative_term2 = 0
    iterative_term3 = 0
    iterative_term4 = 0
    shape_factor = 0
    highly_concave_area = 0
    highly_convex_area = 0

    # Color vertices and calculate metrics
    for i, vertex in enumerate(vertex_tuples):
        # Calculate distance-based factor
        shape_factor += ((get_distance(vertex, inscribed_center)) - inscribed_radius) * (vertex_areas[i]) / inscribed_radius

        # Highly concave regions (red)
        if classifications['mean'][vertex] < 0 and -mean_values[vertex] > (1/inscribed_radius):
            mesh.visual.vertex_colors[i] = [255, 0, 0, 0]
            weighted_concave += -mean_values[vertex] * (vertex_areas[i]) * ((max_distance/2) - (get_distance(vertex, circumsphere.centroid)))
            iterative_term1 += (vertex_areas[i]) * (1 + (1/(mean_values[vertex]*inscribed_radius)))
            iterative_term3 += (vertex_areas[i]) * -mean_values[vertex] * inscribed_radius
            highly_concave_area += vertex_areas[i]

        # Convex regions (yellow)
        if classifications['mean'][vertex] == 1:
            mesh.visual.vertex_colors[i] = [255, 255, 0, 0]

        # Highly convex regions (blue)
        if classifications['mean'][vertex] == 2:
            mesh.visual.vertex_colors[i] = [0, 0, 255, 0]
            highly_convex_area += vertex_areas[i]
            weighted_convex += mean_values[vertex] * (vertex_areas[i]) * (get_distance(vertex, inscribed_center) - inscribed_radius)
            iterative_term2 += (vertex_areas[i]) * (1 - (1/(mean_values[vertex]*inscribed_radius)))
            iterative_term4 += (vertex_areas[i]) * mean_values[vertex] * inscribed_radius

    # Combined terms
    iterative_term = iterative_term1 + iterative_term2
    iterative_term_nl = iterative_term3 + iterative_term4
    total_curved_area = highly_concave_area + highly_convex_area

    return {
        'mesh': mesh,
        'weighted_concave': weighted_concave,
        'weighted_convex': weighted_convex,
        'shape_factor': shape_factor,
        'iterative_term': iterative_term,
        'iterative_term1': iterative_term1,
        'iterative_term2': iterative_term2,
        'iterative_term_nl': iterative_term_nl,
        'highly_concave_area': highly_concave_area,
        'highly_convex_area': highly_convex_area,
        'total_curved_area': total_curved_area,
        'mean_roundness': classifications['roundness']['mean']
    }

def calculate_shape_metrics(mesh, convex_hull, area_stats, visualization_results, inscribed_radius, circumsphere):
    """Calculate comprehensive shape metrics."""
    # Volume metrics
    mesh_volume = mesh.volume
    convex_hull_volume = convex_hull.volume
    circumsphere_volume = circumsphere.volume

    # Area metrics
    total_area = mesh.area
    concave_area = area_stats['concave']
    convex_area = area_stats['convex']
    flat_area = area_stats['flat']

    # Calculate theoretical sphere areas
    max_distance = find_maximum_distance(circumsphere)
    inscribed_sphere_area = 4 * math.pi * (inscribed_radius ** 2)
    circumscribed_sphere_area = 4 * math.pi * ((max_distance/2) ** 2)

    # Get visualization metrics
    weighted_concave = visualization_results['weighted_concave']
    weighted_convex = visualization_results['weighted_convex']
    iterative_term = visualization_results['iterative_term']
    iterative_term1 = visualization_results['iterative_term1']
    iterative_term2 = visualization_results['iterative_term2']
    iterative_term_nl = visualization_results['iterative_term_nl']
    shape_factor = visualization_results['shape_factor']
    highly_concave_area = visualization_results['highly_concave_area']
    highly_convex_area = visualization_results['highly_convex_area']
    total_curved_area = visualization_results['total_curved_area']

    # Ratios and combined metrics
    iterative_term_tc_ratio = iterative_term / total_curved_area if total_curved_area > 0 else 0
    shape_factor_ratio = shape_factor / total_area
    iterative_term_ratio = iterative_term / total_area
    iterative_term_nl_ratio = iterative_term_nl / total_area

    return {
        'mean_roundness': visualization_results['mean_roundness'],
        'weighted_concave': weighted_concave,
        'weighted_convex': weighted_convex,
        'total_curved_area': total_curved_area,
        'weighted_relative': (weighted_concave + weighted_convex) / total_area,
        'highly_concave_area': highly_concave_area,
        'highly_convex_area': highly_convex_area,
        'concave_area': concave_area,
        'concave_ratio': concave_area / total_area,
        'convex_area': convex_area,
        'convex_ratio': convex_area / total_area,
        'flat_area': flat_area,
        'inscribed_radius': inscribed_radius,
        'inscribed_sphere_area': inscribed_sphere_area,
        'circumsphere_volume': circumsphere_volume,
        'circumsphere_area': circumscribed_sphere_area,
        'convex_hull_volume': convex_hull_volume,
        'convex_hull_area': convex_hull.area,
        'mesh_volume': mesh_volume,
        'mesh_area': total_area,
        'iterative_term_tc_ratio': iterative_term_tc_ratio,
        'shape_factor_ratio': shape_factor_ratio,
        'iterative_term_ratio': iterative_term_ratio,
        'iterative_term1': iterative_term1,
        'iterative_term2': iterative_term2,
        'iterative_term_nl_ratio': iterative_term_nl_ratio
    }


def preprocess_mesh(mesh):
    
    # Repair mesh
    trimesh.repair.fill_holes(mesh)

    # Create convex hull
    convex_hull = trimesh.convex.convex_hull(mesh, 'Fx')

    return mesh, convex_hull



def analyze_particle_shape(mesh, visualization=True, save_stl = True):
    """Main function to analyze a 3D particle shape."""
    # Load and preprocess mesh
    mesh, convex_hull = preprocess_mesh(mesh)

    # Get mesh components
    faces = mesh.faces
    vertices = mesh.vertices
    circumsphere = mesh.bounding_sphere

    print("Started Analysis!")

    # Calculate curvatures
    principal_curvatures, vertex_areas = get_curvature_and_areas(faces, vertices)

    # Find maximum inscribed sphere
    inscribed_center, inscribed_radius = find_maximum_inscribed_sphere(vertices)

    classifications = classify_curvature(principal_curvatures, vertices, inscribed_radius)

    area_stats = calculate_area_statistics(vertices, classifications, vertex_areas)

    # Visualize mesh and calculate additional metrics
    visualization_results = visualize_curvature(
        mesh, vertices, classifications,
        classifications['mean_values'], inscribed_radius,
        inscribed_center, vertex_areas, circumsphere
    )

    segmented_mesh, segments = region_segmentation(mesh, vertices, classifications)
    
    metrics = calculate_shape_metrics(
        mesh, convex_hull, area_stats, visualization_results,
        inscribed_radius, circumsphere
    )
    
    # Print summary
    print("=== Shape Analysis Results ===")
    print(f"Mesh Volume: {metrics['mesh_volume']:.6f}")
    print(f"Mesh Area: {metrics['mesh_area']:.6f}")
    print(f"Inscribed Radius: {metrics['inscribed_radius']:.6f}")

    print("Done Analysis")

    return segmented_mesh, visualization_results, metrics
