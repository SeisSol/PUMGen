#include "MPIConvenience.h"
#include <mpi.h>
#include <vector>

void sparseAlltoallv(const void* sendbuf, const int* sendsize, const std::size_t* senddisp,
                     MPI_Datatype sendtype, void* recvbuf, const int* recvsize,
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

  std::vector<MPI_Request> requests(commsize * 2, MPI_REQUEST_NULL);

  // (no special handling of self-to-self comm at the moment)

  for (int i = 0; i < commsize; ++i) {
    if (sendsize[i] > 0) {
      MPI_Isend(reinterpret_cast<const char*>(sendbuf) + senddisp[i] * sendtypesize, sendsize[i],
                sendtype, i, Tag, comm, &requests[i]);
    }
  }
  for (int i = 0; i < commsize; ++i) {
    if (recvsize[i] > 0) {
      MPI_Irecv(reinterpret_cast<char*>(recvbuf) + recvdisp[i] * recvtypesize, recvsize[i],
                recvtype, i, Tag, comm, &requests[i + commsize]);
    }
  }

  MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
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
