// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#ifndef PUMGEN_SRC_HELPER_MPICONVENIENCE_H_
#define PUMGEN_SRC_HELPER_MPICONVENIENCE_H_

#include <cstddef>
#include <mpi.h>

void sparseAlltoallv(const void* sendbuf, const std::size_t* sendsize, const std::size_t* senddisp,
                     MPI_Datatype sendtype, void* recvbuf, const std::size_t* recvsize,
                     const std::size_t* recvdisp, MPI_Datatype recvtype, MPI_Comm comm);

void sparseAlltoallv(const void* sendbuf, const int* sendsize, const std::size_t* senddisp,
                     MPI_Datatype sendtype, void* recvbuf, const int* recvsize,
                     const std::size_t* recvdisp, MPI_Datatype recvtype, MPI_Comm comm);

void sparseAlltoallv(const void* sendbuf, const int* sendsize, const int* senddisp,
                     MPI_Datatype sendtype, void* recvbuf, const int* recvsize, const int* recvdisp,
                     MPI_Datatype recvtype, MPI_Comm comm);

void largeScatterv(const void* sendbuf, const std::size_t* sendsize, const std::size_t* senddisp,
                   MPI_Datatype sendtype, void* recvbuf, std::size_t recvsize,
                   MPI_Datatype recvtype, int root, MPI_Comm comm);

#endif // PUMGEN_SRC_HELPER_MPICONVENIENCE_H_
