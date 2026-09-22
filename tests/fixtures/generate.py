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


def main():
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("General.NumThreads", 1)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)
    generate_layered()
    generate_coarse()
    generate_periodic()
    generate_tiny()
    gmsh.finalize()


if __name__ == "__main__":
    main()
