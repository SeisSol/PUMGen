# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause

"""Generate the mesh fixtures and the reference PUML files for the regression tests.

Run from this directory: python3 generate.py
Requires the gmsh Python API, numpy and h5py. The reference files are computed
from the gmsh model directly, i.e. independently of PUMGen.
"""

import itertools
import os
import struct

import gmsh
import h5py
import numpy as np

# local face numbering of a tetrahedron as used by SeisSol
FACE_NODES = [(1, 0, 2), (0, 1, 3), (1, 2, 3), (2, 0, 3)]

WIDTH, DEPTH, LAYER = 2.0, 1.0, 0.25
# OpenCASCADE enlarges bounding boxes slightly, so classify with a generous tolerance
TOL = 1e-3
FREE_SURFACE, INTERFACE, ABSORBING = 101, 103, 105


def boundary_offset(tag):
    return tag - 100 if tag >= 100 else tag


def build_model(size_top, size_bottom):
    gmsh.model.add("layered")
    top = gmsh.model.occ.addBox(0, 0, -LAYER, WIDTH, WIDTH, LAYER)
    bottom = gmsh.model.occ.addBox(0, 0, -DEPTH, WIDTH, WIDTH, DEPTH - LAYER)
    gmsh.model.occ.fragment([(3, top)], [(3, bottom)])
    gmsh.model.occ.synchronize()

    for dim, tag in gmsh.model.getEntities(3):
        zmin = gmsh.model.getBoundingBox(dim, tag)[2]
        gmsh.model.addPhysicalGroup(3, [tag], 1 if zmin > -LAYER - TOL else 2)

    groups = {FREE_SURFACE: [], INTERFACE: [], ABSORBING: []}
    for dim, tag in gmsh.model.getEntities(2):
        _, _, zmin, _, _, zmax = gmsh.model.getBoundingBox(dim, tag)
        if zmin > -TOL:
            groups[FREE_SURFACE].append(tag)
        elif abs(zmin + LAYER) < TOL and abs(zmax + LAYER) < TOL:
            groups[INTERFACE].append(tag)
        else:
            groups[ABSORBING].append(tag)
    for physical, tags in groups.items():
        gmsh.model.addPhysicalGroup(2, tags, physical)

    for dim, tag in gmsh.model.getEntities(0):
        z = gmsh.model.getValue(dim, tag, [])[2]
        gmsh.model.mesh.setSize([(dim, tag)], size_bottom if z < -LAYER - TOL else size_top)

    gmsh.model.mesh.generate(3)


def physical_of(dim, tag):
    physicals = gmsh.model.getPhysicalGroupsForEntity(dim, tag)
    return int(physicals[0]) if len(physicals) > 0 else 0


def reference_data():
    """Mesh data in the layout PUMGen writes: vertices by node tag, cells in file order."""
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    order = np.argsort(node_tags)
    node_tags = node_tags[order]
    assert np.array_equal(node_tags, np.arange(1, len(node_tags) + 1)), "expected dense node tags"
    vertices = coords.reshape(-1, 3)[order]

    face_bc = {}
    for dim, tag in gmsh.model.getEntities(2):
        physical = physical_of(dim, tag)
        if physical == 0:
            continue
        types, _, nodes = gmsh.model.mesh.getElements(dim, tag)
        for element_type, element_nodes in zip(types, nodes):
            _, _, _, num_nodes, _, _ = gmsh.model.mesh.getElementProperties(element_type)
            for tri in np.asarray(element_nodes).reshape(-1, num_nodes):
                face_bc[frozenset(int(v) for v in tri[:3])] = boundary_offset(physical)

    cells, groups, boundary = [], [], []
    for dim, tag in gmsh.model.getEntities(3):
        physical = physical_of(dim, tag)
        types, _, nodes = gmsh.model.mesh.getElements(dim, tag)
        for element_type, element_nodes in zip(types, nodes):
            _, _, _, num_nodes, _, _ = gmsh.model.mesh.getElementProperties(element_type)
            for tet in np.asarray(element_nodes).reshape(-1, num_nodes):
                code = 0
                for face, local in enumerate(FACE_NODES):
                    bc = face_bc.get(frozenset(int(tet[i]) for i in local), 0)
                    code |= bc << (8 * face)
                cells.append(tet - 1)
                groups.append(physical)
                boundary.append(code)
    return vertices, np.asarray(cells), np.asarray(groups), np.asarray(boundary)


def periodic_identification(num_vertices):
    """Each vertex is identified with the smallest vertex index of its periodic class."""
    parent = list(range(num_vertices))

    def find(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v

    for dim in range(3):
        for _, tag in gmsh.model.getEntities(dim):
            master, nodes, master_nodes, _ = gmsh.model.mesh.getPeriodicNodes(dim, tag)
            if master == tag:
                continue
            for node, master_node in zip(nodes, master_nodes):
                a, b = find(int(node) - 1), find(int(master_node) - 1)
                parent[max(a, b)] = min(a, b)
    return np.array([find(v) for v in range(num_vertices)])


def write_reference(path, mesh_file, periodic=False):
    """Reference PUML data from a written mesh file, so coordinates match their text form."""
    gmsh.clear()
    gmsh.open(mesh_file)
    vertices, cells, groups, boundary = reference_data()
    with h5py.File(path, "w") as out:
        out.create_dataset("connect", data=cells.astype("<u8"))
        out.create_dataset("geometry", data=vertices.astype("<f8"))
        out.create_dataset("group", data=groups.astype("<i4"))
        out.create_dataset("boundary", data=boundary.astype("<i4"))
        if periodic:
            out.create_dataset("identify", data=periodic_identification(len(vertices)).astype("<u8"))


def read_msh41_ascii(path):
    """The nodes by tag and the node tags of the three-dimensional elements in the order of the
    file, read from an ASCII MSH 4.1 file without gmsh."""
    lines = open(path).read().split("\n")
    nodes = {}
    cells = []
    i = 0
    while i < len(lines):
        if lines[i] == "$Nodes":
            blocks = int(lines[i + 1].split()[0])
            i += 2
            for _ in range(blocks):
                count = int(lines[i].split()[3])
                tags = [int(lines[i + 1 + k]) for k in range(count)]
                for k in range(count):
                    nodes[tags[k]] = [float(x) for x in lines[i + 1 + count + k].split()[:3]]
                i += 1 + 2 * count
        elif lines[i] == "$Elements":
            blocks = int(lines[i + 1].split()[0])
            i += 2
            for _ in range(blocks):
                dim, _, _, count = (int(x) for x in lines[i].split())
                if dim == 3:
                    cells.extend([int(x) for x in lines[i + 1 + k].split()[1:]] for k in range(count))
                i += 1 + count
        else:
            i += 1
    return nodes, cells


def write_high_order_reference(path, mesh_file, linear_reference):
    """Reference of a mesh of order 2 whose vertices are the ones of the linear reference: the
    other nodes of every cell, in the order of gmsh, with their offsets and the orders."""
    nodes, cells = read_msh41_ascii(mesh_file)
    with h5py.File(linear_reference, "r") as linear:
        connect = linear["connect"][...]
        datasets = {name: linear[name][...] for name in ("connect", "geometry", "group", "boundary")}
    # the vertices are the nodes which are vertices of a cell, numbered in the order of their tags;
    # they have to give the linear reference
    vertex = {tag: i for i, tag in enumerate(sorted({tag for cell in cells for tag in cell[:4]}))}
    assert all([vertex[tag] for tag in cell[:4]] == list(connect[c]) for c, cell in enumerate(cells))
    geometry = [x for cell in cells for tag in cell[4:] for x in nodes[tag]]
    offsets = np.cumsum([0] + [3 * (len(cell) - 4) for cell in cells])
    with h5py.File(path, "w") as out:
        for name, data in datasets.items():
            out.create_dataset(name, data=data)
        out.create_dataset("geometry_ho", data=np.array(geometry, dtype="<f8"))
        out["geometry_ho"].attrs.create("node-ordering", "gmsh", dtype=h5py.string_dtype("ascii"))
        out.create_dataset("geometry_ho_offsets", data=offsets.astype("<u8"))
        out.create_dataset("order", data=np.full(len(cells), 2, dtype="u1"))


def write_mesh(path, version, binary=False):
    gmsh.option.setNumber("Mesh.MshFileVersion", version)
    gmsh.option.setNumber("Mesh.Binary", 1 if binary else 0)
    gmsh.write(path)
    gmsh.option.setNumber("Mesh.Binary", 0)


# nodes per element type, for the element types in the fixtures
ELEMENT_NODES = {15: 1, 1: 2, 2: 3, 4: 4, 8: 3, 9: 6, 11: 10}


def convert_binary(source, target, byte_order, size_bytes):
    """Rewrites a binary MSH 4.1 file (as written by gmsh: little endian, 8-byte sizes) with
    another byte order ("<" or ">") and size of size_t (4 or 8)."""
    with open(source, "rb") as f:
        data = f.read()
    pos = 0
    out = bytearray()
    size_format = "Q" if size_bytes == 8 else "I"

    def line():
        nonlocal pos
        end = data.index(b"\n", pos)
        text = data[pos:end]
        pos = end + 1
        return text

    def get(fmt, count=1):
        nonlocal pos
        fmt = f"<{count}{fmt}"
        values = struct.unpack_from(fmt, data, pos)
        pos += struct.calcsize(fmt)
        out.extend(struct.pack(byte_order + fmt[1:].replace("Q", size_format), *values))
        return values

    def sizes(count=1):
        return get("Q", count)

    assert line() == b"$MeshFormat" and line() == b"4.1 1 8"
    out.extend(f"$MeshFormat\n4.1 1 {size_bytes}\n".encode())
    assert get("i")[0] == 1
    assert line() == b""
    out.extend(b"\n")
    while pos < len(data):
        text = line()
        out.extend(text + b"\n")
        if text == b"$MeshFormat":
            continue
        if text == b"$Entities":
            counts = sizes(4)
            for kind, count in enumerate(counts):
                for _ in range(count):
                    get("i")
                    get("d", 3 if kind == 0 else 6)
                    get("i", sizes()[0])
                    if kind > 0:
                        get("i", sizes()[0])
        elif text == b"$Nodes":
            num_blocks = sizes(4)[0]
            for _ in range(num_blocks):
                dim, _, parametric = get("i", 3)
                count = sizes()[0]
                sizes(count)
                get("d", count * (3 + (dim if parametric else 0)))
        elif text == b"$Elements":
            num_blocks = sizes(4)[0]
            for _ in range(num_blocks):
                _, _, element_type = get("i", 3)
                count = sizes()[0]
                sizes(count * (1 + ELEMENT_NODES[element_type]))
        elif text == b"$Periodic":
            for _ in range(sizes()[0]):
                get("i", 3)
                get("d", sizes()[0])
                sizes(2 * sizes()[0])
        elif text.startswith(b"$End"):
            continue
        else:
            raise ValueError(f"unexpected section {text}")
        assert line() == b""
        end = line()
        assert end == b"$End" + text[1:]
        out.extend(b"\n" + end + b"\n")
    with open(target, "wb") as f:
        f.write(out)


def counts_exercise_chunking(count):
    # the element/vertex counts are chosen such that the default chunk distribution differs
    # from a ceil-based one for 3 and for 4 ranks
    return count % 3 == 1 and count % 4 in (1, 2)


def generate_layered():
    for size_top in itertools.count(0.2, 0.0025):
        gmsh.clear()
        build_model(size_top, 2.0 * size_top)
        vertices, cells, _, _ = reference_data()
        if counts_exercise_chunking(len(cells)) and counts_exercise_chunking(len(vertices)):
            break
    print(f"layered: {len(cells)} cells, {len(vertices)} vertices, size {size_top:.4f}")

    write_mesh("layered-v41.msh", 4.1)
    write_mesh("layered-binary-v41.msh", 4.1, binary=True)
    write_mesh("layered-v22.msh", 2.2)
    gmsh.write("layered.neu")
    write_reference("layered.puml.h5", "layered-v41.msh")


def generate_coarse():
    gmsh.clear()
    build_model(0.5, 1.0)
    vertices, cells, _, _ = reference_data()
    print(f"coarse: {len(cells)} cells, {len(vertices)} vertices")

    write_mesh("coarse-v41.msh", 4.1)
    write_mesh("coarse-binary-v41.msh", 4.1, binary=True)

    # dense node tags that do not start at 1
    node_tags, _, _ = gmsh.model.mesh.getNodes()
    gmsh.model.mesh.renumberNodes(node_tags, node_tags + 1000)
    write_mesh("coarse-offset-v41.msh", 4.1)
    gmsh.model.mesh.renumberNodes(node_tags + 1000, node_tags)

    gmsh.model.mesh.setOrder(2)
    write_mesh("coarse-o2-v41.msh", 4.1)
    write_mesh("coarse-o2-binary-v41.msh", 4.1, binary=True)
    write_mesh("coarse-o2-v22.msh", 2.2)
    gmsh.model.mesh.setOrder(1)

    gmsh.model.mesh.partition(2)
    write_mesh("coarse-partitioned-v41.msh", 4.1)

    # the converter reproduces gmsh's layout exactly before it is used for other layouts
    convert_binary("coarse-binary-v41.msh", "coarse-binary-copy.msh", "<", 8)
    with open("coarse-binary-v41.msh", "rb") as a, open("coarse-binary-copy.msh", "rb") as b:
        assert a.read() == b.read()
    os.remove("coarse-binary-copy.msh")
    convert_binary("coarse-binary-v41.msh", "coarse-binary-bigendian-v41.msh", ">", 8)
    convert_binary("coarse-binary-v41.msh", "coarse-binary-size4-v41.msh", "<", 4)
    with open("coarse-binary-v41.msh", "rb") as f:
        data = f.read()
    with open("coarse-binary-truncated-v41.msh", "wb") as f:
        f.write(data[: data.index(b"$Elements") + 400])

    write_reference("coarse.puml.h5", "coarse-v41.msh")
    write_reference("coarse-o2.puml.h5", "coarse-o2-v41.msh")


def generate_periodic():
    gmsh.clear()
    gmsh.model.add("periodic")
    gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
    gmsh.model.occ.synchronize()
    translation = [1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    left = [tag for _, tag in gmsh.model.getEntitiesInBoundingBox(-TOL, -TOL, -TOL, TOL, 1 + TOL, 1 + TOL, 2)]
    right = [tag for _, tag in gmsh.model.getEntitiesInBoundingBox(1 - TOL, -TOL, -TOL, 1 + TOL, 1 + TOL, 1 + TOL, 2)]
    gmsh.model.mesh.setPeriodic(2, right, left, translation)
    gmsh.model.addPhysicalGroup(3, [1], 1)
    others = [tag for _, tag in gmsh.model.getEntities(2) if tag not in left + right]
    gmsh.model.addPhysicalGroup(2, others, FREE_SURFACE)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.4)
    gmsh.model.mesh.generate(3)
    write_mesh("periodic-v41.msh", 4.1)
    write_mesh("periodic-binary-v41.msh", 4.1, binary=True)
    gmsh.option.setNumber("Mesh.SaveParametric", 1)
    write_mesh("periodic-parametric-v41.msh", 4.1)
    gmsh.option.setNumber("Mesh.SaveParametric", 0)
    write_reference("periodic.puml.h5", "periodic-v41.msh", periodic=True)


def generate_periodic_tiny():
    """A periodic cube of as few vertices as gmsh makes, to be converted on more ranks than it has
    vertices."""
    gmsh.clear()
    gmsh.model.add("periodic-tiny")
    gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
    gmsh.model.occ.synchronize()
    translation = [1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    left = [tag for _, tag in gmsh.model.getEntitiesInBoundingBox(-TOL, -TOL, -TOL, TOL, 1 + TOL, 1 + TOL, 2)]
    right = [tag for _, tag in gmsh.model.getEntitiesInBoundingBox(1 - TOL, -TOL, -TOL, 1 + TOL, 1 + TOL, 1 + TOL, 2)]
    gmsh.model.mesh.setPeriodic(2, right, left, translation)
    gmsh.model.addPhysicalGroup(3, [1], 1)
    others = [tag for _, tag in gmsh.model.getEntities(2) if tag not in left + right]
    gmsh.model.addPhysicalGroup(2, others, FREE_SURFACE)
    # the eight corners only
    for _, curve in gmsh.model.getEntities(1):
        gmsh.model.mesh.setTransfiniteCurve(curve, 2)
    for _, surface in gmsh.model.getEntities(2):
        gmsh.model.mesh.setTransfiniteSurface(surface)
    gmsh.model.mesh.setTransfiniteVolume(1)
    gmsh.model.mesh.generate(3)
    write_mesh("periodic-tiny-v41.msh", 4.1)
    write_mesh("periodic-tiny-binary-v41.msh", 4.1, binary=True)
    write_reference("periodic-tiny.puml.h5", "periodic-tiny-v41.msh", periodic=True)


TINY_V41 = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Entities
0 0 0 1
1 0 0 0 1 1 1 1 7 0
$EndEntities
$Nodes
1 4 1 4
3 1 0 4
1
2
3
4
0 0 0
1 0 0
0 1 0
0 0 1
$EndNodes
$Elements
1 1 1 1
3 1 4 1
1 1 2 3 4
$EndElements
"""

TINY_V22 = """$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
4
1 0 0 0
2 1 0 0
3 0 1 0
4 0 0 1
$EndNodes
$Elements
1
1 4 2 7 1 1 2 3 4
$EndElements
"""


def tiny_binary(node_tags, element_nodes, node_header=(1, 4, 1, 4)):
    """The single tetrahedron of TINY_V41 as binary MSH 4.1 file."""
    data = b"$MeshFormat\n4.1 1 8\n" + struct.pack("<i", 1) + b"\n$EndMeshFormat\n"
    data += b"$Entities\n" + struct.pack("<4Q", 0, 0, 0, 1)
    data += struct.pack("<i6dQiQ", 1, 0, 0, 0, 1, 1, 1, 1, 7, 0) + b"\n$EndEntities\n"
    data += b"$Nodes\n" + struct.pack("<4Q", *node_header) + struct.pack("<3iQ", 3, 1, 0, 4)
    data += struct.pack("<4Q", *node_tags)
    data += struct.pack("<12d", 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1) + b"\n$EndNodes\n"
    data += b"$Elements\n" + struct.pack("<4Q", 1, 1, 1, 1) + struct.pack("<3iQ", 3, 1, 4, 1)
    data += struct.pack("<5Q", 1, *element_nodes) + b"\n$EndElements\n"
    return data


def generate_tiny():
    variants = {
        "tiny-v41.msh": TINY_V41,
        "tiny-nophysical-v41.msh": TINY_V41.replace("1 0 0 0 1 1 1 1 7 0", "1 0 0 0 1 1 1 0 0"),
        "tiny-v40.msh": TINY_V41.replace("4.1 0 8", "4 0 8"),
        "tiny-truncated-v41.msh": TINY_V41.split("$Elements")[0] + "$Elements\n1 1 1 1\n3 1 4 1\n1 1 2\n",
        "tiny-badnumber-v41.msh": TINY_V41.replace("1 0 0\n0 1 0", "1 abc 0\n0 1 0"),
        "tiny-unknowntype-v41.msh": TINY_V41.replace("3 1 4 1\n", "3 1 999 1\n"),
        "tiny-sparse-v41.msh": TINY_V41.replace("1 4 1 4\n", "1 4 1 5\n")
        .replace("3\n4\n0 0 0", "3\n5\n0 0 0")
        .replace("1 1 2 3 4\n", "1 1 2 3 5\n"),
        "tiny-duplicate-v41.msh": TINY_V41.replace("1\n2\n3\n4\n", "1\n2\n2\n4\n"),
        "tiny-badnode-v41.msh": TINY_V41.replace("1 1 2 3 4\n", "1 1 2 3 9\n"),
        "tiny-v22.msh": TINY_V22,
        "tiny-duplicate-v22.msh": TINY_V22.replace("3 0 1 0", "2 0 1 0"),
        "tiny-badnode-v22.msh": TINY_V22.replace("1 2 3 4\n$EndElements", "1 2 3 9\n$EndElements"),
        "tiny-duplicateface-v22.msh": TINY_V22.replace(
            "1\n1 4 2 7 1 1 2 3 4\n", "3\n1 4 2 7 1 1 2 3 4\n2 2 2 101 1 1 3 2\n3 2 2 105 1 2 1 3\n"
        ),
    }
    for name, content in variants.items():
        with open(name, "w") as out:
            out.write(content)
    binary_variants = {
        "tiny-binary-v41.msh": tiny_binary((1, 2, 3, 4), (1, 2, 3, 4)),
        "tiny-binary-duplicate-v41.msh": tiny_binary((1, 2, 2, 4), (1, 2, 3, 4)),
        "tiny-binary-badnode-v41.msh": tiny_binary((1, 2, 3, 4), (1, 2, 3, 9)),
        "tiny-binary-sparse-v41.msh": tiny_binary((1, 2, 3, 5), (1, 2, 3, 5), (1, 4, 1, 5)),
    }
    for name, content in binary_variants.items():
        with open(name, "wb") as out:
            out.write(content)
    write_reference("tiny.puml.h5", "tiny-v41.msh")


# the faces of the cells in the numbering of PUML, as the specification of the file format has them
FACES = {
    "tet": [(1, 0, 2), (0, 1, 3), (1, 2, 3), (2, 0, 3)],
    "hex": [(0, 4, 7, 3), (1, 2, 6, 5), (0, 1, 5, 4), (3, 7, 6, 2), (0, 3, 2, 1), (4, 5, 6, 7)],
    "wedge": [(0, 2, 1), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)],
    "pyramid": [(0, 3, 2, 1), (0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)],
}
MSH_TYPES = {"tet": 4, "hex": 5, "wedge": 6, "pyramid": 7}
FACE_TYPES = {3: 2, 4: 3}
VTK_TYPES = {"tet": 10, "hex": 12, "wedge": 13, "pyramid": 14}
VOLUME_NAMES = {1: "hexahedron", 2: "pyramids", 3: "wedges", 4: "tetrahedra"}
SURFACE_NAMES = {1: "free surface", 3: "front", 5: "absorbing"}


def mixed_cells():
    """A conforming mesh of all four kinds of cells, cube by cube along x: a hexahedron, six
    pyramids around the centre of the next cube, two wedges in the third, and six tetrahedra in a
    cube next to the triangles of the wedges. Returns the coordinates of the nodes and the cells as
    (kind, nodes, group), with the nodes in the order of gmsh and positively oriented."""
    points = []
    index = {}

    def node(x, y, z):
        if (x, y, z) not in index:
            index[(x, y, z)] = len(points)
            points.append((x, y, z))
        return index[(x, y, z)]

    def cube(x, y, z):
        return [node(x + dx, y + dy, z + dz) for dz in (0, 1) for (dx, dy) in ((0, 0), (1, 0), (1, 1), (0, 1))]

    cells = [("hex", cube(0, 0, 0), 1)]
    corners = cube(1, 0, 0)
    apex = node(1.5, 0.5, 0.5)
    for face in FACES["hex"]:
        # the base of a pyramid turns towards its apex, the hexahedral faces point outwards
        cells.append(("pyramid", [corners[v] for v in reversed(face)] + [apex], 2))
    cells.append(("wedge", [node(2, 0, 0), node(2, 0, 1), node(3, 0, 0), node(2, 1, 0), node(2, 1, 1), node(3, 1, 0)], 3))
    cells.append(("wedge", [node(3, 0, 0), node(2, 0, 1), node(3, 0, 1), node(3, 1, 0), node(2, 1, 1), node(3, 1, 1)], 3))
    # Kuhn's six tetrahedra around the diagonal from (3, 1, 0), which splits the face y = 1 as the
    # triangles of the wedges do
    steps = {"u": (-1, 0, 0), "v": (0, 1, 0), "w": (0, 0, 1)}
    for order in itertools.permutations("uvw"):
        position = [3, 1, 0]
        tet = [node(*position)]
        for axis in order:
            position = [a + b for a, b in zip(position, steps[axis])]
            tet.append(node(*position))
        cells.append(("tet", tet, 4))

    for kind, nodes, _ in cells:
        if kind == "tet":
            a, b, c, d = (np.array(points[n]) for n in nodes)
            if np.dot(np.cross(b - a, c - a), d - a) < 0:
                nodes[2], nodes[3] = nodes[3], nodes[2]
    return points, cells


def boundary_faces(cells):
    """The faces of the cells found once, with their boundary condition: 1 on top (z = 1), 3 at
    y = 0, and 5 elsewhere."""
    count = {}
    for kind, nodes, _ in cells:
        for face in FACES[kind]:
            key = tuple(sorted(nodes[v] for v in face))
            count[key] = count.get(key, 0) + 1
    return {key for key, n in count.items() if n == 1}


def bc_of(face, points):
    if all(points[v][2] == 1 for v in face):
        return 1
    if all(points[v][1] == 0 for v in face):
        return 3
    return 5


def generate_mixed():
    points, cells = mixed_cells()
    boundary = boundary_faces(cells)
    gmsh.clear()
    gmsh.model.add("mixed")
    for group in (1, 2, 3, 4):
        gmsh.model.addDiscreteEntity(3, group)
        gmsh.model.addPhysicalGroup(3, [group], group, name=VOLUME_NAMES[group])
    for bc in (1, 3, 5):
        gmsh.model.addDiscreteEntity(2, bc)
        gmsh.model.addPhysicalGroup(2, [bc], 100 + bc, name=SURFACE_NAMES[bc])
    gmsh.model.mesh.addNodes(3, 1, list(range(1, len(points) + 1)), [x for p in points for x in p])
    tag = 1
    for kind, nodes, group in cells:
        gmsh.model.mesh.addElementsByType(group, MSH_TYPES[kind], [tag], [n + 1 for n in nodes])
        tag += 1
    for face in sorted(boundary):
        # the vertices of a face in the cyclic order of the cell that has it
        cyclic = next([nodes[v] for v in f] for kind, nodes, _ in cells for f in FACES[kind]
                      if tuple(sorted(nodes[v] for v in f)) == face)
        gmsh.model.mesh.addElementsByType(bc_of(face, points), FACE_TYPES[len(face)], [tag], [n + 1 for n in cyclic])
        tag += 1
    write_mesh("mixed-v41.msh", 4.1)
    write_mesh("mixed-binary-v41.msh", 4.1, binary=True)
    write_mesh("mixed-v22.msh", 2.2)

    # the references, with the cells in the order of the files: gmsh writes the MSH 2.2 file by
    # element type
    _, order = read_msh41_ascii("mixed-v41.msh")
    write_mixed_reference("mixed.puml.h5", points, cells, boundary, order, read_physical_names("mixed-v41.msh"))
    write_mixed_reference(
        "mixed-v22.puml.h5", points, cells, boundary, read_msh22_cells("mixed-v22.msh"), read_physical_names("mixed-v22.msh")
    )


def read_physical_names(path):
    """The names of the physical groups of an ASCII MSH file, in its order, as (dimension, tag, name)."""
    lines = open(path).read().split("\n")
    first = lines.index("$PhysicalNames") + 2
    names = []
    for line in lines[first : first + int(lines[first - 1])]:
        dimension, tag, name = line.split(" ", 2)
        names.append((int(dimension), int(tag), name.strip().strip('"')))
    return names


def add_names(dataset, names, dimension, boundary):
    """The attributes naming the values of a dataset: ids, as the dataset holds them, and names."""
    chosen = [(tag, name) for dim, tag, name in names if dim == dimension]
    if chosen:
        ids = [tag - 100 if boundary and tag >= 100 else tag for tag, _ in chosen]
        dataset.attrs.create("ids", np.array(ids, dtype="<i4"))
        dataset.attrs.create("names", [name for _, name in chosen], dtype=h5py.string_dtype("utf-8"))


def read_msh22_cells(path):
    """The node tags of the three-dimensional elements of an ASCII MSH 2.2 file, in its order."""
    lines = open(path).read().split("\n")
    first = lines.index("$Elements") + 2
    cells = []
    for line in lines[first : first + int(lines[first - 1])]:
        values = [int(x) for x in line.split()]
        if values[1] in MSH_TYPES.values():
            cells.append(values[3 + values[2] :])
    return cells


def write_mixed_reference(path, points, cells, boundary, order, names):
    by_nodes = {tuple(n + 1 for n in nodes): (kind, nodes, group) for kind, nodes, group in cells}
    ordered = [by_nodes[tuple(tags)] for tags in order]
    faces = np.zeros((len(ordered), 6), dtype="<i4")
    for c, (kind, nodes, _) in enumerate(ordered):
        for f, face in enumerate(FACES[kind]):
            key = tuple(sorted(nodes[v] for v in face))
            if key in boundary:
                faces[c, f] = bc_of(key, points)
    with h5py.File(path, "w") as out:
        out.create_dataset("connect", data=np.array([n for _, nodes, _ in ordered for n in nodes], dtype="<u8"))
        out.create_dataset("connect_offsets", data=np.cumsum([0] + [len(nodes) for _, nodes, _ in ordered]).astype("<u8"))
        out.create_dataset("cell_type", data=np.array([VTK_TYPES[kind] for kind, _, _ in ordered], dtype="u1"))
        out.create_dataset("geometry", data=np.array(points, dtype="<f8"))
        out.create_dataset("group", data=np.array([group for _, _, group in ordered], dtype="<i4"))
        out.create_dataset("boundary", data=faces)
        add_names(out["group"], names, 3, False)
        add_names(out["boundary"], names, 2, True)


def main():
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("General.NumThreads", 1)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)
    generate_layered()
    generate_coarse()
    generate_periodic()
    generate_periodic_tiny()
    generate_tiny()
    generate_mixed()
    gmsh.finalize()


if __name__ == "__main__":
    main()
