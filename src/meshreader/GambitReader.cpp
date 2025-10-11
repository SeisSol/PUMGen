// SPDX-FileCopyrightText: 2017 Technical University of Munich
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileContributor: Sebastian Rettenberger <sebastian.rettenberger@tum.de>

#include "GambitReader.h"

const char* puml::GambitReader::GAMBIT_FILE_ID = "** GAMBIT NEUTRAL FILE";
const char* puml::GambitReader::ENDSECTION = "ENDOFSECTION";
const char* puml::GambitReader::NODAL_COORDINATES = "NODAL COORDINATES";
const char* puml::GambitReader::ELEMENT_CELLS = "ELEMENTS/CELLS";
const char* puml::GambitReader::ELEMENT_GROUP = "ELEMENT GROUP";
const char* puml::GambitReader::BOUNDARY_CONDITIONS = "BOUNDARY CONDITIONS";
