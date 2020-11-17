/**
 * @file
 *  This file is part of PUMGen
 *
 *  For conditions of distribution and use, please see the copyright
 *  notice in the file 'COPYING' at the root directory of this package
 *  and the copyright notice at https://github.com/SeisSol/PUMGen
 *
 * @copyright 2017 Technical University of Munich
 * @author Sebastian Rettenberger <sebastian.rettenberger@tum.de>
 */

#include "GambitReader.h"

const char* puml::GambitReader::GAMBIT_FILE_ID = "** GAMBIT NEUTRAL FILE";
const char* puml::GambitReader::ENDSECTION = "ENDOFSECTION";
const char* puml::GambitReader::NODAL_COORDINATES = "NODAL COORDINATES";
const char* puml::GambitReader::ELEMENT_CELLS = "ELEMENTS/CELLS";
const char* puml::GambitReader::ELEMENT_GROUP = "ELEMENT GROUP";
const char* puml::GambitReader::BOUNDARY_CONDITIONS = "BOUNDARY CONDITIONS";
