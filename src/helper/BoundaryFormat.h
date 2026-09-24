// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#ifndef PUMGEN_SRC_HELPER_BOUNDARYFORMAT_H_
#define PUMGEN_SRC_HELPER_BOUNDARYFORMAT_H_

#include <cstddef>
#include <cstdint>
#include <vector>

#include "mesh/CellType.h"

/**
 * Packs the boundary conditions of the first four faces of every cell into one integer, with
 * bitsPerFace bits per face, the first face in the lowest bits. The faces hold as many conditions
 * per cell as its kind has faces. A condition which does not fit is an error.
 */
std::vector<std::uint64_t> packBoundaries(const std::vector<int>& faces,
                                          const std::vector<puml::CellType>& cellTypes,
                                          int bitsPerFace);

/**
 * The boundary conditions of every cell in facesPerCell columns, 0 behind the faces of a cell.
 */
std::vector<std::int32_t> boundariesPerFace(const std::vector<int>& faces,
                                            const std::vector<puml::CellType>& cellTypes,
                                            std::size_t facesPerCell);

#endif // PUMGEN_SRC_HELPER_BOUNDARYFORMAT_H_
