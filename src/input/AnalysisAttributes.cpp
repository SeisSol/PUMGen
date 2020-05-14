/**
 * @file
 *  This file is part of PUMGen
 *
 *  For conditions of distribution and use, please see the copyright
 *  notice in the file 'COPYING' at the root directory of this package
 *  and the copyright notice at https://github.com/SeisSol/PUMGen
 *
 * @copyright 2018 LMU-TUM 
 * @author Thomas Ulrich
 */

#include "AnalysisAttributes.h"

AnalysisAttributes::AnalysisAttributes(const char* xmlFilename, int numFaces) {
    readXmlFile(xmlFilename);
    faceBound.reserve(numFaces);
}

void AnalysisAttributes::readXmlFile(const char* xmlFilename) {
    if(doc.LoadFile(xmlFilename) !=XML_SUCCESS) {
       logError() << "Unable to open or to parse file" << xmlFilename;
    }
}

void AnalysisAttributes::set_BoundaryConditions() {
    std::string line,sval;
    XMLElement* pRoot;
    XMLElement* child;
    pRoot = doc.FirstChildElement("freeSurface");
    if (pRoot) {
       std::vector<std::string> tokens;
       line =  pRoot-> GetText();
       split(tokens, line, ',');
       for(int i = 0; i < tokens.size(); i++) {
           faceBound.emplace_back( faceBoundary(std::atoi(tokens[i].c_str())-1, 1));
       }
    }
    pRoot = doc.FirstChildElement("dynamicRupture");
    if (pRoot) {
       std::vector<std::string> tokens;
       line =  pRoot-> GetText();
       split(tokens, line, ',');
       for(int i = 0; i < tokens.size(); i++) {
           faceBound.emplace_back( faceBoundary(std::atoi(tokens[i].c_str())-1, 3));
       }
    }
    pRoot = doc.FirstChildElement("absorbing");
    if (pRoot) {
       std::vector<std::string> tokens;
       line =  pRoot-> GetText();
       split(tokens, line, ',');
       for(int i = 0; i < tokens.size(); i++) {
           faceBound.emplace_back( faceBoundary(std::atoi(tokens[i].c_str())-1, 5));
       }
    }
    //more general way of defining boundaryCondition
    child = doc.FirstChildElement("boundaryCondition");
    if (child) {
       for (child; child; child = child->NextSiblingElement("boundaryCondition"))
       {
          sval = child->Attribute("tag");
          int faultTag =  atoi(sval.c_str());
          std::vector<std::string> tokens;
          line =  child-> GetText();
          split(tokens, line, ',');
          for(int i = 0; i < tokens.size(); i++) {
              faceBound.emplace_back( faceBoundary(std::atoi(tokens[i].c_str())-1, faultTag));
          }
       }
    }
}
