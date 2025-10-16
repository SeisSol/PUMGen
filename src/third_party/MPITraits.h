// SPDX-FileCopyrightText: 2025 SeisSol Group
// SPDX-FileCopyrightText: 2020 Ludwig-Maximilians-Universität München
//
// SPDX-License-Identifier: BSD-3-Clause
#ifndef PUMGEN_SRC_THIRD_PARTY_MPITRAITS_H_
#define PUMGEN_SRC_THIRD_PARTY_MPITRAITS_H_

#include <mpi.h>

namespace tndm {

template <typename T> struct mpi_type;
template <typename T> inline MPI_Datatype mpi_type_t() { return mpi_type<T>::type(); }

// integer types
template <> struct mpi_type<int> {
  static MPI_Datatype type() { return MPI_INT; }
};
template <> struct mpi_type<long> {
  static MPI_Datatype type() { return MPI_LONG; }
};
template <> struct mpi_type<long long> {
  static MPI_Datatype type() { return MPI_LONG_LONG; }
};
template <> struct mpi_type<unsigned> {
  static MPI_Datatype type() { return MPI_UNSIGNED; }
};
template <> struct mpi_type<unsigned long> {
  static MPI_Datatype type() { return MPI_UNSIGNED_LONG; }
};
template <> struct mpi_type<unsigned long long> {
  static MPI_Datatype type() { return MPI_UNSIGNED_LONG_LONG; }
};

// floating point types
template <> struct mpi_type<float> {
  static MPI_Datatype type() { return MPI_FLOAT; }
};
template <> struct mpi_type<double> {
  static MPI_Datatype type() { return MPI_DOUBLE; }
};

template <typename T> class mpi_array_type {
  public:
  mpi_array_type(int count) {
    MPI_Type_contiguous(count, mpi_type_t<T>(), &type);
    MPI_Type_commit(&type);
  }

  ~mpi_array_type() { MPI_Type_free(&type); }

  MPI_Datatype const& get() const { return type; }

  private:
  MPI_Datatype type;
};

} // namespace tndm

#endif // PUMGEN_SRC_THIRD_PARTY_MPITRAITS_H_
