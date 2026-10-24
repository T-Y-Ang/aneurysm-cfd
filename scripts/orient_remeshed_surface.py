# -------------------------------------------------------------
# CORRECT ORIENTATION OF REMESHED ANEURYSM SURFACE
# -------------------------------------------------------------
#
# Project: Patient-specific intracranial aneurysm CFD
#
# Purpose:
# Correct inconsistent triangle orientations in the remeshed
# aneurysm surface before attempting tetrahedral volume meshing.
#
# The original patient-derived geometry will remain unchanged.
# -------------------------------------------------------------

from pathlib import Path
from collections import defaultdict, deque

import gmsh
import numpy as np
import struct


# -------------------------------------------------------------
# DEFINE INPUT AND OUTPUT PATHS
# -------------------------------------------------------------

# Input: surface mesh generated during our remeshing experiment.
input_path = Path(
    "mesh/patient_aneurysm_surface_remesh_test.msh"
)

# Output: a separate STL containing corrected triangle orientations.
output_path = Path(
    "geometry/patient_aneurysm_remeshed_oriented.stl"
)


# -------------------------------------------------------------
# INITIALISE GMSH AND LOAD THE SURFACE
# -------------------------------------------------------------

gmsh.initialize()

try:
    gmsh.open(str(input_path))

    print("\nREMESHED ANEURYSM SURFACE")
    print("Input file:", input_path)

    # Retrieve all node identifiers and their XYZ coordinates.
    # The original mesh coordinates are expressed in millimetres.
    node_tags, coordinates, _ = gmsh.model.mesh.getNodes()

    points = np.asarray(
        coordinates,
        dtype=np.float64
    ).reshape(-1, 3)

    # Create a lookup table:
    # Gmsh node identifier -> XYZ coordinates.
    #
    # This is necessary because Gmsh node identifiers are not
    # guaranteed to match their positions in the points array.
    node_lookup = {
        int(tag): points[i]
        for i, tag in enumerate(node_tags)
    }

    print("Number of nodes:", len(node_tags))
    print("Coordinate units: millimetres")

        # -------------------------------------------------------------
    # EXTRACT TRIANGULAR SURFACE ELEMENTS
    # -------------------------------------------------------------
    #
    # Gmsh meshes may contain different element types:
    # - Points (0D)
    # - Lines (1D)
    # - Triangles (2D)
    # - Tetrahedra (3D)
    #
    # We only need the 2D triangular elements because we are
    # correcting the orientation of the aneurysm surface.
    # -------------------------------------------------------------

    element_types, _, element_nodes = (
        gmsh.model.mesh.getElements(dim=2)
    )

    triangle_blocks = []

    for element_type, connectivity in zip(
        element_types,
        element_nodes
    ):

        # Gmsh element type 2 represents a three-node triangle.
        if element_type == 2:

            triangles_in_block = np.asarray(
                connectivity,
                dtype=np.int64
            ).reshape(-1, 3)

            triangle_blocks.append(triangles_in_block)

    # Combine triangles from all Gmsh surface entities.
    if not triangle_blocks:
        raise RuntimeError(
            "No three-node surface triangles were found."
        )

    triangles = np.vstack(triangle_blocks)

    print("\nSURFACE TRIANGLE EXTRACTION")
    print("Number of triangles:", len(triangles))
    print("Triangle array shape:", triangles.shape)
    print("First triangle node IDs:", triangles[0])

        # -------------------------------------------------------------
    # IDENTIFY SHARED TRIANGLE EDGES
    # -------------------------------------------------------------
    #
    # Each triangle has three directed edges:
    # A -> B, B -> C, and C -> A.
    #
    # For a correctly oriented surface, neighbouring triangles
    # traverse their shared edge in opposite directions.
    # -------------------------------------------------------------

    edge_owners = defaultdict(list)

    for triangle_index, (a, b, c) in enumerate(triangles):

        for u, v in [(a, b), (b, c), (c, a)]:

            u = int(u)
            v = int(v)

            # Sort node IDs to identify the same physical edge,
            # regardless of which direction it is traversed.
            edge_key = tuple(sorted((u, v)))

            # Record the direction of the edge.
            direction = 1 if (u, v) == edge_key else -1

            edge_owners[edge_key].append(
                (triangle_index, direction)
            )

    # Check whether every edge belongs to exactly two triangles.
    open_edges = sum(
        len(owners) == 1
        for owners in edge_owners.values()
    )

    nonmanifold_edges = sum(
        len(owners) > 2
        for owners in edge_owners.values()
    )

    # Count edges whose two triangles traverse them
    # in the same direction.
    inconsistent_edges = sum(
        owners[0][1] == owners[1][1]
        for owners in edge_owners.values()
        if len(owners) == 2
    )

    print("\nSHARED EDGE ANALYSIS")
    print("Unique edges:", len(edge_owners))
    print("Open edges:", open_edges)
    print("Non-manifold edges:", nonmanifold_edges)
    print("Inconsistently oriented edges:", inconsistent_edges)

    # Stop if the surface topology is unsuitable for correction.
    if open_edges != 0 or nonmanifold_edges != 0:
        raise RuntimeError(
            "Surface is not a closed two-manifold mesh."
        )

        # -------------------------------------------------------------
    # BUILD TRIANGLE ADJACENCY RELATIONSHIPS
    # -------------------------------------------------------------
    #
    # Two triangles are neighbours if they share an edge.
    #
    # For every neighbouring pair, determine whether one
    # triangle must be reversed relative to the other.
    # -------------------------------------------------------------

    adjacency = [[] for _ in range(len(triangles))]

    for owners in edge_owners.values():

        (triangle_i, direction_i), (
            triangle_j, direction_j
        ) = owners

        # If both triangles traverse their shared edge in the
        # same direction, exactly one must be reversed.
        must_flip = direction_i == direction_j

        adjacency[triangle_i].append(
            (triangle_j, must_flip)
        )

        adjacency[triangle_j].append(
            (triangle_i, must_flip)
        )

    # -------------------------------------------------------------
    # PROPAGATE ORIENTATION CORRECTIONS
    # -------------------------------------------------------------
    #
    # Breadth-first search (BFS) visits neighbouring triangles
    # and assigns each one a required orientation.
    #
    # False = retain original vertex order.
    # True  = reverse vertex order.
    # -------------------------------------------------------------

    flips = {}
    contradictions = 0
    connected_components = 0

    for start in range(len(triangles)):

        if start in flips:
            continue

        connected_components += 1

        # Begin each connected component with an unchanged
        # reference triangle.
        flips[start] = False

        queue = deque([start])

        while queue:

            current = queue.popleft()

            for neighbour, must_flip in adjacency[current]:

                # XOR determines whether the neighbour must
                # have the opposite flip state.
                required_flip = (
                    flips[current] ^ must_flip
                )

                if neighbour not in flips:

                    flips[neighbour] = required_flip
                    queue.append(neighbour)

                elif flips[neighbour] != required_flip:

                    # Conflicting requirements indicate that
                    # consistent orientation was not possible.
                    contradictions += 1

    # -------------------------------------------------------------
    # REPORT RESULTS
    # -------------------------------------------------------------

    triangles_to_flip = sum(flips.values())

    print("\nORIENTATION CORRECTION PLAN")
    print("Connected components:", connected_components)
    print("Triangles requiring reversal:", triangles_to_flip)
    print("Orientation contradictions:", contradictions)

    if connected_components != 1:
        raise RuntimeError(
            "Expected one connected aneurysm surface."
        )

    if contradictions != 0:
        raise RuntimeError(
            "Surface orientation constraints are contradictory."
        )

    print("PASS: Consistent triangle orientation is possible.")

        # -------------------------------------------------------------
    # APPLY TRIANGLE ORIENTATION CORRECTIONS
    # -------------------------------------------------------------
    #
    # Work on a copy so the original triangle array remains
    # unchanged for comparison and debugging.
    # -------------------------------------------------------------

    corrected_triangles = triangles.copy()

    for triangle_index, should_flip in flips.items():

        if should_flip:

            # Reversing two vertices reverses the triangle normal.
            corrected_triangles[
                triangle_index, [1, 2]
            ] = corrected_triangles[
                triangle_index, [2, 1]
            ]

    print("\nORIENTATION CORRECTIONS APPLIED")
    print("Triangles reversed:", triangles_to_flip)


    # -------------------------------------------------------------
    # CALCULATE ENCLOSED VOLUME
    # -------------------------------------------------------------
    #
    # Each oriented surface triangle forms a signed tetrahedron
    # with a common reference point.
    #
    # Summing these signed contributions gives the enclosed
    # volume for a consistently oriented closed surface.
    # -------------------------------------------------------------

    reference_point = points.mean(axis=0)

    signed_volume = 0.0

    for a, b, c in corrected_triangles:

        p0 = node_lookup[int(a)] - reference_point
        p1 = node_lookup[int(b)] - reference_point
        p2 = node_lookup[int(c)] - reference_point

        contribution = np.dot(
            p0,
            np.cross(p1, p2)
        ) / 6.0

        signed_volume += contribution

    # -------------------------------------------------------------
    # ENSURE THE SURFACE NORMALS POINT OUTWARDS
    # -------------------------------------------------------------
    #
    # For a consistently oriented closed surface:
    # Positive signed volume indicates outward orientation.
    # Negative signed volume indicates inward orientation.
    # -------------------------------------------------------------

    if signed_volume < 0:

        corrected_triangles[:, [1, 2]] = (
            corrected_triangles[:, [2, 1]]
        )

        signed_volume = -signed_volume

        print("Reversed entire surface to face outwards.")

    # -------------------------------------------------------------
    # COMPARE AGAINST ORIGINAL PATIENT-DERIVED GEOMETRY
    # -------------------------------------------------------------

    original_volume = 170.073934

    volume_difference_percent = (
        (signed_volume - original_volume)
        / original_volume
        * 100
    )

    print("\nCORRECTED SURFACE VOLUME CHECK")
    print("Original volume (mm³):", original_volume)
    print("Corrected volume (mm³):", signed_volume)
    print(
        "Volume difference (%):",
        volume_difference_percent
    )

    # A provisional volume-preservation threshold.
    # This is a diagnostic criterion, not proof of CFD quality.
    if abs(volume_difference_percent) > 2.0:

        raise RuntimeError(
            "Corrected surface differs from original "
            "by more than 2% in enclosed volume."
        )

    print("PASS: Enclosed volume is within 2% of original.")

        # -------------------------------------------------------------
    # EXPORT CORRECTED SURFACE AS BINARY STL
    # -------------------------------------------------------------
    #
    # Binary STL stores each triangle using:
    # - One surface normal vector
    # - Three XYZ vertex coordinates
    # - One two-byte attribute field
    #
    # Coordinates remain in millimetres.
    # -------------------------------------------------------------

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    with open(output_path, "wb") as stl_file:

        # STL header: exactly 80 bytes.
        header = (
            b"Patient aneurysm - orientation corrected"
        )

        stl_file.write(header.ljust(80, b"\0"))

        # Number of triangular facets.
        stl_file.write(
            struct.pack("<I", len(corrected_triangles))
        )

        for a, b, c in corrected_triangles:

            p0 = node_lookup[int(a)]
            p1 = node_lookup[int(b)]
            p2 = node_lookup[int(c)]

            # Calculate the outward-facing triangle normal.
            normal = np.cross(
                p1 - p0,
                p2 - p0
            )

            normal_length = np.linalg.norm(normal)

            if normal_length == 0:
                raise RuntimeError(
                    "Degenerate triangle encountered during export."
                )

            normal = normal / normal_length

            # STL facet format:
            # 12 float32 values followed by one uint16.
            facet_data = [
                *normal,
                *p0,
                *p1,
                *p2
            ]

            stl_file.write(
                struct.pack("<12fH", *facet_data, 0)
            )

    print("\nCORRECTED STL EXPORT")
    print("Output file:", output_path)
    print("Triangles exported:", len(corrected_triangles))
    print(
        "File size (bytes):",
        output_path.stat().st_size
    )
    print("PASS: Corrected STL exported successfully.")



finally:
    # Always release Gmsh resources, even if loading fails.
    gmsh.finalize()