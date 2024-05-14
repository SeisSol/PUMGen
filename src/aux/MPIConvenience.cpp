#include "MPIConvenience.h"
#include <limits>
#include <mpi.h>
#include <vector>

namespace {
constexpr std::size_t DataIncrement = std::numeric_limits<int>::max();
constexpr int Tag = 1000;

std::vector<MPI_Request> largeIsend(const void* sendbuf, std::size_t sendsize, std::size_t senddisp,
                                    MPI_Datatype sendtype, int dest, int tag, MPI_Comm comm) {
  int sendtypesizePre;
  std::vector<MPI_Request> requests;
  MPI_Type_size(sendtype, &sendtypesizePre);
  std::size_t sendtypesize = sendtypesizePre;
  for (std::size_t position = 0; position < sendsize; position += DataIncrement) {
    int localSize = static_cast<int>(std::min(sendsize - position, DataIncrement));
    requests.push_back(MPI_REQUEST_NULL);
    MPI_Isend(reinterpret_cast<const char*>(sendbuf) + (senddisp + position) * sendtypesize,
              localSize, sendtype, dest, tag, comm, &requests.back());
  }
  return requests;
}

std::vector<MPI_Request> largeIrecv(void* recvbuf, std::size_t recvsize, std::size_t recvdisp,
                                    MPI_Datatype recvtype, int dest, int tag, MPI_Comm comm) {
  int recvtypesizePre;
  std::vector<MPI_Request> requests;
  MPI_Type_size(recvtype, &recvtypesizePre);
  std::size_t recvtypesize = recvtypesizePre;
  for (std::size_t position = 0; position < recvsize; position += DataIncrement) {
    int localSize = static_cast<int>(std::min(recvsize - position, DataIncrement));
    requests.push_back(MPI_REQUEST_NULL);
    MPI_Irecv(reinterpret_cast<char*>(recvbuf) + (recvdisp + position) * recvtypesize, localSize,
              recvtype, dest, tag, comm, &requests.back());
  }
  return requests;
}
} // namespace

void sparseAlltoallv(const void* sendbuf, const std::size_t* sendsize, const std::size_t* senddisp,
                     MPI_Datatype sendtype, void* recvbuf, const std::size_t* recvsize,
                     const std::size_t* recvdisp, MPI_Datatype recvtype, MPI_Comm comm) {
  int commsize;
  int commrank;
  int sendtypesizePre;
  int recvtypesizePre;
  MPI_Comm_size(comm, &commsize);
  MPI_Comm_rank(comm, &commrank);

  std::vector<MPI_Request> requests;
  requests.reserve(commsize * 2);

  // (no special handling of self-to-self comm at the moment)

  for (int i = 0; i < commsize; ++i) {
    auto localRequests = largeIsend(sendbuf, sendsize[i], senddisp[i], sendtype, i, Tag, comm);
    requests.insert(requests.end(), localRequests.begin(), localRequests.end());
  }
  for (int i = 0; i < commsize; ++i) {
    auto localRequests = largeIrecv(recvbuf, recvsize[i], recvdisp[i], recvtype, i, Tag, comm);
    requests.insert(requests.end(), localRequests.begin(), localRequests.end());
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

void largeScatterv(const void* sendbuf, const std::size_t* sendsize, const std::size_t* senddisp,
                   MPI_Datatype sendtype, void* recvbuf, std::size_t recvsize,
                   MPI_Datatype recvtype, int root, MPI_Comm comm) {
  std::vector<MPI_Request> requests;
  int commrank;
  int commsize;
  MPI_Comm_size(comm, &commsize);
  MPI_Comm_rank(comm, &commrank);
  if (commrank == root) {
    requests.reserve(commsize + 1);
    for (int i = 0; i < commsize; ++i) {
      auto sendRequests = largeIsend(sendbuf, sendsize[i], senddisp[i], sendtype, i, Tag, comm);
      requests.insert(requests.end(), sendRequests.begin(), sendRequests.end());
    }
  }
  auto recvRequests = largeIrecv(recvbuf, recvsize, 0, recvtype, root, Tag, comm);
  requests.insert(requests.end(), recvRequests.begin(), recvRequests.end());
  MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
}
