#ifndef PUMGEN_EASIMESHSIZE_H
#define PUMGEN_EASIMESHSIZE_H

#include <array>
#include <memory>
#include <string>
#include "easi/Component.h"

// Workaround for easi linking bug
namespace easi {
class YAMLParser;
}
struct ElasticMaterial {
  double lambda, mu, rho;
};

class EasiMeshSize {
  private:
  easi::YAMLParser* parser;
  easi::Component* model; // Unique ptr to model leads to segfault
  std::string easiFileName;
  public:
  EasiMeshSize();;
  EasiMeshSize(std::string easiFileName);

  double getMeshSize(std::array<double, 3> point);

  ~EasiMeshSize();
};

#endif // PUMGEN_EASIMESHSIZE_H
