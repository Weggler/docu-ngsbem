from math import dist

import pandas as pd


def max_mesh_size(mesh):
    """Return the maximum edge length of the surface mesh."""
    points = mesh.ngmesh.Points()
    h = 0.0

    for element in mesh.ngmesh.Elements2D():
        vertices = [vertex.nr for vertex in element.vertices]
        for first, second in zip(vertices, vertices[1:] + vertices[:1]):
            h = max(h, dist(tuple(points[first]), tuple(points[second])))

    return h


def append_results(results, csv_path):
    write_header = not csv_path.exists()
    pd.DataFrame(results).to_csv(
        csv_path,
        mode="w" if write_header else "a",
        header=write_header,
        index=False,
    )
