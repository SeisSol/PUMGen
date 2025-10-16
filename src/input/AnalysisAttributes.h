// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#ifndef PUMGEN_SRC_INPUT_ANALYSISATTRIBUTES_H_
#define PUMGEN_SRC_INPUT_ANALYSISATTRIBUTES_H_

#include "split.h"
#include "tinyxml2/tinyxml2.h"
#include "utils/logger.h"
#include <list>
#include <string>
#include <unordered_map>
#include <vector>

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

#endif // PUMGEN_SRC_INPUT_ANALYSISATTRIBUTES_H_
