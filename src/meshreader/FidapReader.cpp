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

#include "FidapReader.h"

const char* FidapReader::FIDAP_FILE_ID = "** FIDAP NEUTRAL FILE";
const char* FidapReader::NODAL_COORDINATES = "NODAL COORDINATES";
const char* FidapReader::BOUNDARY_CONDITIONS = "BOUNDARY CONDITIONS";
const char* FidapReader::ELEMENT_GROUPS = "ELEMENT GROUPS";

const char* FidapReader::ZONE_GROUP = "GROUP:";
const char* FidapReader::ZONE_ELEMENTS = "ELEMENTS:";
const char* FidapReader::ZONE_NODES = "NODES:";
const char* FidapReader::ZONE_GEOMETRY = "GEOMETRY:";
const char* FidapReader::ZONE_TYPE = "TYPE:";
