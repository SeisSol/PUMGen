#ifndef PUMGEN_EASIMESHSIZE_H
#define PUMGEN_EASIMESHSIZE_H

#include <array>
#include <string_view>

class EasiMeshSize {
  public:
  EasiMeshSize(std::string_view easiFileName) {
    // TODO(Lukas) Do something
  }

  double getMeshSize(std::array<double, 3> point) {
    if(point[2] > 0 ) {
      return 1.0;
    } else {
      return 0.1;
    }
  }
};

#endif // PUMGEN_EASIMESHSIZE_H
