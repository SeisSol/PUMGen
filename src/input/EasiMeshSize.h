#ifndef PUMGEN_EASIMESHSIZE_H
#define PUMGEN_EASIMESHSIZE_H

#include "easi/Component.h"
#include <MeshTypes.h>
#include <SimModel.h>
#include <array>
#include <memory>
#include <string>

// Workaround for easi linking bug
namespace easi {
class YAMLParser;
}
struct ElasticMaterial {
  double lambda, mu, rho;
};

class EasiMeshSize {
private:
  std::string easiFileName;
  double targetedFrequency;
  double elementsPerWaveLength;
  easi::YAMLParser* parser;
  easi::Component* model; // Unique ptr to model leads to segfault
  pGModel simModel;

  int findGroup(std::array<double, 3> point) {
    GRIter regionIt = GM_regionIter(simModel);
    while (pGRegion region = GRIter_next(regionIt)) {
      // Note: Does not work on discrete regions!
      if (GR_containsPoint(region, point.data()) > 0) {
        return GEN_tag(region);
      }

      return -1;
    }
  }

  public:
  EasiMeshSize();;
  EasiMeshSize(std::string easiFileName, double targetedFrequency, double elementsPerWaveLength,
               pGModel simModel);

  double getMeshSize(std::array<double, 3> point);

  ~EasiMeshSize();
};

#endif // PUMGEN_EASIMESHSIZE_H
