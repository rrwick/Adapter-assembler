// Copyright 2017 Ryan Wick

// This file is part of Adapter-assembler

// Adapter-assembler is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
// version.

// Adapter-assembler is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.

// You should have received a copy of the GNU General Public License along with Adapter-assembler.  If not, see
// <http://www.gnu.org/licenses/>.


#include "misc.h"

#include <iostream>
#include <sstream>
#include <iomanip>


std::string int_to_string(long long n) {
    std::stringstream ss;
    ss.imbue(std::locale(""));
    ss << std::fixed << n;
    return ss.str();
}


void print_hash_progress(std::string filename, long long base_count) {
    std::cerr << "\r  " << filename << " (" << int_to_string(base_count) << " bp)";
}
