// This is a -*- C++ -*- header file.

/* print.h -- Print data structures

   Copyright 2002-2004 Jesus A. De Loera, David Haws, Raymond
      Hemmecke, Peter Huggins, Jeremy Tauzer, Ruriko Yoshida
   Copyright 2006 Matthias Koeppe

   This file is part of LattE.
   
   LattE is free software; you can redistribute it and/or modify it
   under the terms of the version 2 of the GNU General Public License
   as published by the Free Software Foundation.

   LattE is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with LattE; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
*/

#ifndef PRINT_H
#define PRINT_H

#include "cone.h"
#include "cone_consumer.h"

void printVector(const vec_ZZ &, int);
void printListVector(listVector*, int);
listVector* Changeproj(listVector* basis, listVector* equ, int numOfVars);
void printRationalVector(rationalVector*, int);
void printCone(listCone*, int);
void printListCone(listCone*, int);
void printVectorToFile(ostream &, const vec_ZZ &, int);
void printListVectorToFile(ostream &, listVector*, int);
void printVectorToFileWithoutBrackets(ostream &, const vec_ZZ &, int);
void printListVectorToFileWithoutBrackets(ostream &, listVector*, int);
void printRationalVectorToFile(ostream &, rationalVector*, int);
void printRationalVectorToFileWithoutBrackets(ostream &, rationalVector*, 
					      int);
void printConeToFile(ostream &out, listCone* cones, int numOfVars);

/* ConeFileHeader is declared in cone.h, so that cone_consumer.h can use it
   too.  Reading, writing and checking it live here. */
void printConeFileHeaderToFile(ostream &out, const ConeFileHeader &header);
ConeFileHeader readConeFileHeader(istream &in);
void checkConeFileHeader(const ConeFileHeader &header, const ConeFileHeader &expected,
		    const string &filename, const string &used_option);

listCone *
readConeFromFile(istream &in);

void printListConeToFile(const char*, listCone*, int);

// Likewise, prefixed by a header block describing the cones.
void printListConeToFile(const char*, listCone*, int, const ConeFileHeader &);

// Read a list of cones in the format of `printListCone'.
// If HEADER_OUT is non-NULL and the file has a header block, it is stored there.
listCone *
readListConeFromStream(istream &in, ConeFileHeader *header_out = NULL);

listCone *
readListConeFromFile(const char *, ConeFileHeader *header_out = NULL);

// Likewise, but feed the cones one by one to CONSUMER.
void
readListConeFromFile(istream &in, ConeConsumer &consumer,
		     ConeFileHeader *header_out = NULL);

void printResidueFile(const char*, listCone*, int);

void
print_debug_matrix(const mat_ZZ &);
void
print_debug_vector(const vec_ZZ &);

#endif
