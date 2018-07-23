
#include <algorithm>
#include <cstring>
#include <string>

#include <list>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>


std::vector<std::string> &split(std::vector<std::string> &elems, const std::string &s, char delim) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


