#ifndef ANALYSISATTRIBUTES_H
#define ANALYSISATTRIBUTES_H

#include <unordered_map>
#include "tinyxml2/tinyxml2.h"
#include <string>
#include <list>
#include <vector>
#include "utils/logger.h"
#include "split.h"



using namespace tinyxml2;

struct faceBoundary {
  faceBoundary(int faceID, int bcType) : faceID(faceID), bcType(bcType) {};
  int faceID;
  int bcType;
};

class AnalysisAttributes {
  public:
    std::vector<faceBoundary> faceBound;
    AnalysisAttributes(const char* xmlFilename, int numFaces);
  private:
    XMLDocument doc;
    void readXmlFile(const char* xmlFilename);
    void set_BoundaryConditions();
};

#endif
