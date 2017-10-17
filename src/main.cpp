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


#include <iostream>
#include <vector>
#include <unordered_map>

#include "arguments.h"
#include "kmers.h"

#define PROGRAM_VERSION "0.1.0"


int main(int argc, char **argv)
{
    Arguments args(argc, argv);
    if (args.parsing_result == BAD)
        return 1;
    else if (args.parsing_result == HELP)
        return 0;
    else if (args.parsing_result == VERSION) {
        std::cout << "Adapter-assembler v" << PROGRAM_VERSION << "\n";
        return 0;
    }

    std::cerr << "\n";

    Kmers kmers(args.kmer);
    kmers.add_fastq(args.input_reads, args.start, args.margin);
    kmers.remove_low_depth_kmers(args.min_depth);
    kmers.output_gfa();

    std::cerr << "\n";
    return 0;
}
