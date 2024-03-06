#include "MPIConvenience.h"
#include <limits>
#include <mpi.h>
#include <vector>

void sparseAlltoallv(const void* sendbuf, const std::size_t* sendsize, const std::size_t* senddisp,
                     MPI_Datatype sendtype, void* recvbuf, const std::size_t* recvsize,
                     const std::size_t* recvdisp, MPI_Datatype recvtype, MPI_Comm comm) {
  int commsize;
  int commrank;
  int sendtypesizePre;
  int recvtypesizePre;
  MPI_Comm_size(comm, &commsize);
  MPI_Comm_rank(comm, &commrank);
  MPI_Type_size(sendtype, &sendtypesizePre);
  MPI_Type_size(sendtype, &recvtypesizePre);

  constexpr int Tag = 1000;

  std::size_t sendtypesize = sendtypesizePre;
  std::size_t recvtypesize = recvtypesizePre;

  std::vector<MPI_Request> requests;
  requests.reserve(commsize * 2);

  // (no special handling of self-to-self comm at the moment)

  constexpr std::size_t DataIncrement = std::numeric_limits<int>::max();

  for (int i = 0; i < commsize; ++i) {
    for (std::size_t position = 0; position < sendsize[i]; position += DataIncrement) {
      int localSize = static_cast<int>(std::min(sendsize[i] - position, DataIncrement));
      requests.push_back(MPI_REQUEST_NULL);
      MPI_Isend(reinterpret_cast<const char*>(sendbuf) + senddisp[i] * sendtypesize + position,
                localSize, sendtype, i, Tag, comm, &requests.back());
    }
  }
  for (int i = 0; i < commsize; ++i) {
    for (std::size_t position = 0; position < recvsize[i]; position += DataIncrement) {
      int localSize = static_cast<int>(std::min(recvsize[i] - position, DataIncrement));
      requests.push_back(MPI_REQUEST_NULL);
      MPI_Irecv(reinterpret_cast<char*>(recvbuf) + recvdisp[i] * recvtypesize + position, localSize,
                recvtype, i, Tag, comm, &requests.back());
    }
  }

  MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
}

void sparseAlltoallv(const void* sendbuf, const int* sendsize, const std::size_t* senddisp,
                     MPI_Datatype sendtype, void* recvbuf, const int* recvsize,
                     const std::size_t* recvdisp, MPI_Datatype recvtype, MPI_Comm comm) {
  int commsize;
  MPI_Comm_size(comm, &commsize);
  std::vector<std::size_t> sendsizeSize(sendsize, sendsize + commsize);
  std::vector<std::size_t> recvsizeSize(recvsize, recvsize + commsize);
  sparseAlltoallv(sendbuf, sendsizeSize.data(), senddisp, sendtype, recvbuf, recvsizeSize.data(),
                  recvdisp, recvtype, comm);
}

void sparseAlltoallv(const void* sendbuf, const int* sendsize, const int* senddisp,
                     MPI_Datatype sendtype, void* recvbuf, const int* recvsize, const int* recvdisp,
                     MPI_Datatype recvtype, MPI_Comm comm) {
  int commsize;
  MPI_Comm_size(comm, &commsize);
  std::vector<std::size_t> senddispSize(senddisp, senddisp + commsize);
  std::vector<std::size_t> recvdispSize(recvdisp, recvdisp + commsize);
  sparseAlltoallv(sendbuf, sendsize, senddispSize.data(), sendtype, recvbuf, recvsize,
                  recvdispSize.data(), recvtype, comm);
}
