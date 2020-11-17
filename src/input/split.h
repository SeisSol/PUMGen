#ifndef INPUT_SPLIT_H
#define INPUT_SPLIT_H

#include <algorithm>
#include <cstring>
#include <string>

#include <fstream>
#include <iostream>
#include <list>
#include <vector>

std::vector<std::string> &split(std::vector<std::string> &elems,
                                const std::string &s, char delim);

#endif
