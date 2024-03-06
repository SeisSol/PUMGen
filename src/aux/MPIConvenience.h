#ifndef PUMGEN_AUX_MPI_CONVENIENCE_H_
#define PUMGEN_AUX_MPI_CONVENIENCE_H_

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

#endif // PUMGEN_AUX_MPI_CONVENIENCE_H_
