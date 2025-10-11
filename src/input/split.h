// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#ifndef PUMGEN_SRC_INPUT_SPLIT_H_
#define PUMGEN_SRC_INPUT_SPLIT_H_

#include <algorithm>
#include <cstring>
#include <string>

#include <fstream>
#include <iostream>
#include <list>
#include <vector>

std::vector<std::string>& split(std::vector<std::string>& elems, const std::string& s, char delim);

#endif // PUMGEN_SRC_INPUT_SPLIT_H_
