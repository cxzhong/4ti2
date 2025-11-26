/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2006 4ti2 team.
Main author(s): Raymond Hemmecke.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA. 
*/

void printVersionInfo();
void printVector(vector, int);
void printListVector(listVector*, int);
void printVectorToFile(FILE*, vector, int);
void printListVectorToFile(const char*, listVector*, int);
void printSubsetOfListVectorToFile(const char*, listVector*, vector, int);
void printListVectorWithoutColumnToFile(const char*, listVector*, int, int);
void printListVectorWithAdditionalColumnToFile(const char*, listVector*, int,
					       int, int);
void printTransposedListVectorToFile(const char*, listVector*, int);
void printVectorToFileWithBrackets(FILE*, vector, int);
void printBinomialToFile(FILE*, vector, int, char**);
void printListBinomialsToFile(const char*, listVector*, int, char**);
void printMonomialToFile(FILE*, vector, int, char**);
void printListMonomialsAndBinomialsToFile(const char*, listVector*, int,
					  char**);
void printMatrix(vector, int, int);
void printVectorToFileMaple(FILE*, vector, int);
void printListVectorMaple(const char*, listVector*, int);
void printListVectorMacaulay2(const char*, listVector*, int);
void print3wayTables(const char*, listVector*, int, int, int, int);
void printL1NormOfListVector(listVector*, int);
void printListVectorWithGivenDegreesToFile(const char*, listVector*, int,
					   int, int);
void printSupportsOfListVector(listVector*, int);
void printTypesOfListVector(listVector*, int, int);
void printListVectorWithGivenSupportsToFile(const char*, listVector*, int,
					    int, int);
void printListVectorWithGivenNonzeroEntryToFile(const char*, listVector*, int,
						int);
void writeResult(listVector*, int, const char*, const char*, int);
void printListRepresentativesToFile(const char*, listOrbit*, int);
void printRationalVector(rationalVector*, int);
void printRationalVectorToFileWithoutBrackets(FILE*, rationalVector*, int);

