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

#ifndef ARGUMENTS_H
#define ARGUMENTS_H


#include <string>
#include <vector>


enum ParsingResult {GOOD, BAD, HELP, VERSION};


class Arguments
{
public:
    Arguments(int argc, char **argv);

    ParsingResult parsing_result;

    std::string input_reads;

    int kmer;
    double filter_depth;
    int margin;
    bool start;
    bool end;


private:
    bool does_file_exist(std::string fileName);
};

#endif // ARGUMENTS_H
