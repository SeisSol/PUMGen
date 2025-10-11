// SPDX-FileCopyrightText: 2017 Technical University of Munich
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileContributor: Sebastian Rettenberger <sebastian.rettenberger@tum.de>

#include <mpi.h>

#include "ParallelVertexFilter.h"

MPI_Datatype ParallelVertexFilter::vertexType = MPI_DATATYPE_NULL;
