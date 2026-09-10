/* print.cpp -- Print data structures

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

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <time.h>
#include <list>
#include <vector>

#include "cone.h"
#include "ramon.h"
#include "print.h"
#include "LattException.h"
#include "gnulib/pathmax.h"

/* ----------------------------------------------------------------- */
void printVector(const vec_ZZ &v, int numOfVars)
{
  printVectorToFile(cout, v, numOfVars);
}
/* ----------------------------------------------------------------- */
void printListVector(listVector* basis, int numOfVars)
{
  printListVectorToFile(cout, basis, numOfVars);
}
/* ----------------------------------------------------------------- */
void printRationalVector(rationalVector *v, int numOfVars) {
  int i;

//    if (v==0) {
//      cout << "[]\n";
//      return ;
//    }

  cout << "[";
  for (i=0; i<(numOfVars-1); i++) {
    if ((v->denominators())[i]==1)
      cout << (v->numerators())[i] << " ";
    else
      cout << (v->numerators())[i] << "/" << (v->denominators())[i] << " ";
  }

  if ((v->denominators())[i]==1)
    cout << (v->numerators())[i] << "]" << endl;
  else
    cout << (v->numerators())[i] << "/" << (v->denominators())[i] << "]" << endl;
  return ;
}
/* ----------------------------------------------------------------- */
void printCone(listCone* cones, int numOfVars) {
  printConeToFile(cout, cones, numOfVars);
}
/* ----------------------------------------------------------------- */
void printListCone(listCone* cones, int numOfVars) {
  if (cones==0) cout << "No cones in list.\n";
  while(cones) {
    printCone(cones,numOfVars);
    cones = cones->rest;
  }
  cout << endl;
  return ;
}
/* ----------------------------------------------------------------- */
void printVectorToFile(ostream & out, const vec_ZZ &v, int numOfVars)
{
  int i;
  assert(v.length() == numOfVars);
  out << "[";
  for (i=0; i<(numOfVars-1); i++) {
    out << v[i] << " ";
  }
  out << v[i] << "]\n";
  return ;
}
/* ----------------------------------------------------------------- */
void printListVectorToFile(ostream & out, listVector* basis, int numOfVars) {
  if (basis==0) {
    out << "[]\n";
    return;
  }

  while(basis) {
    printVectorToFile(out,basis->first,numOfVars);
    basis = basis->rest;
  }
  return ;
}
/* ----------------------------------------------------------------- */
void printVectorToFileWithoutBrackets(ostream & out, const vec_ZZ &v, 
				      int numOfVars) {
  int i;

//    if (v==0) return ;

  for (i=0; i<(numOfVars-1); i++) {
    out << v[i] << " ";
  }
  out << v[i] << endl;
  return ;
}
/* ----------------------------------------------------------------- */
void printListVectorToFileWithoutBrackets(ostream & out, listVector* basis, 
					  int numOfVars) {
  if (basis==0) {
    out << numOfVars << " 0\n";
    return;
  }

  while(basis) {
    printVectorToFileWithoutBrackets(out,basis->first,numOfVars);
    basis = basis->rest;
  }
  return ;
}
/* ----------------------------------------------------------------- */
void printRationalVectorToFile(ostream & out, rationalVector *v, 
			       int numOfVars) {
  int i;

  if (v==0) {
    out << "[]\n";
    return ;
  }
  out << "[";
  for (i=0; i<(numOfVars-1); i++) {
    if ((v->denominators())[i]==1)
      out << (v->numerators())[i] << " ";
    else
      out << (v->numerators())[i] << "/" << (v->denominators())[i] << " ";
  }
  
  if ((v->denominators())[i]==1)
    out << (v->numerators())[i] << "]\n";
  else
    out << (v->numerators())[i] << "/" << (v->denominators())[i] << "]\n";
  return ;
}
/* ----------------------------------------------------------------- */
void printRationalVectorToFileWithoutBrackets(ostream & out, 
					      rationalVector *v, 
					      int numOfVars) {
  int i;

  if (v==0) {
    return ;
  }
  for (i=0; i<(numOfVars); i++) {
    if ((v->denominators())[i]==1)
      out << (v->numerators())[i] << " "; 
    else
      out << (v->numerators())[i] << "/" << (v->denominators())[i] << " ";
  }

  out << endl;
  return ;
}
/* ----------------------------------------------------------------- */
void printConeToFile(ostream & out,listCone* cones, int numOfVars)
{
  out << "==========\n";
  out << "Cone.\n";

  out << "Coefficient: " << cones->coefficient << endl;

  out << "Vertex: ";
  printRationalVectorToFile(out,cones->vertex->vertex,numOfVars);

  out << "Extreme rays:\n";  
  printListVectorToFile(out,cones->rays,numOfVars);

  out << "Determinant:" << cones->determinant << endl;
  
  out << "Facets:\n";  
  printListVectorToFile(out,cones->facets,numOfVars);

  out << "Dual determinant:" << cones->dual_determinant << endl;
  
  out << "Lattice points in parallelepiped:\n";
  printListVectorToFile(out,cones->latticePoints,numOfVars);
  out << "==========\n\n";

  return ;
}
/* ----------------------------------------------------------------- */


static bool
look_for(istream &in, const char *token)
{
  string s;
  while (in.good()) {
    in >> s;
    if (s == token) return true;
  }
  return false;
}

/* Scan forward to TOKEN, picking up a `Determinant:VALUE' field on the way and
   storing it in DET. DET is left untouched when no determinant is seen. */
static bool
look_for_capturing_determinant(istream &in, const char *token, ZZ &det)
{
  string s;
  const string prefix = "Determinant:";
  while (in.good()) {
    in >> s;
    if (s == token) return true;
    if (s.compare(0, prefix.size(), prefix) == 0) {
      if (s.size() > prefix.size()) {
	istringstream value(s.substr(prefix.size()));
	value >> det;
      }
      else
	in >> det;
    }
  }
  return false;
}

static void
skip_space(istream &in)
{
  while (isspace(in.peek())) {
    char c;
    in.get(c);
  }
}

/* ----------------------------------------------------------------- */
/* The header block.  See the format description below. */

void
printConeFileHeaderToFile(ostream &out, const ConeFileHeader &header)
{
  out << "LattE cone file.\n";
  out << "Format version: " << header.version << endl;
  out << "Ambient dimension: " << header.numOfVars << endl;
  out << "Homogenized: " << (header.homogenized ? 1 : 0) << endl;
  out << "Dualized: " << (header.dualized ? 1 : 0) << endl;
  out << "Unbounded: " << (header.unbounded ? 1 : 0) << endl;
}

static string
trim(const string &s)
{
  size_t b = s.find_first_not_of(" \t\r\n");
  if (b == string::npos) return "";
  size_t e = s.find_last_not_of(" \t\r\n");
  return s.substr(b, e - b + 1);
}

static bool
parse_header_bool(const string &value, bool &result)
{
  if (value == "0" || value == "false") { result = false; return true; }
  if (value == "1" || value == "true") { result = true; return true; }
  return false;
}

static void
header_parse_error(const string &what)
{
  cerr << "Malformed header in cone file: " << what << "." << endl;
  THROW_LATTE(LattException::fe_Parse, 0);
}

/* How a header describes itself in a diagnostic. */
static string
describe_cones(const ConeFileHeader &header)
{
  if (!header.homogenized) return "vertex cones";
  if (header.dualized) return "dual homogenized cones";
  return "primal homogenized cones";
}

ConeFileHeader
readConeFileHeader(istream &in)
{
  ConeFileHeader header;
  header.present = false;
  skip_space(in);
  /* A file without a header starts with `==========' or `No cones in list.',
     never with the `L' of `LattE cone file.'. */
  if (in.peek() != 'L') return header;

  istream::pos_type start = in.tellg();
  string line;
  if (!getline(in, line)) return header;
  if (trim(line) != "LattE cone file.") {
    if (start == istream::pos_type(-1))
      header_parse_error("the file begins with `" + trim(line)
			 + "' where a cone or a header was expected");
    in.clear();
    in.seekg(start);
    return header;
  }

  header.present = true;
  bool seen_version = false;
  for (;;) {
    istream::pos_type line_start = in.tellg();
    if (!getline(in, line)) break;
    string field = trim(line);
    if (field.size() == 0) continue;
    size_t colon = field.find(':');
    if (colon == string::npos) {
      /* `==========' or `No cones in list.' -- the header is over. */
      in.clear();
      in.seekg(line_start);
      break;
    }
    string key = trim(field.substr(0, colon));
    string value = trim(field.substr(colon + 1));
    if (key == "Format version") {
      istringstream s(value);
      if (!(s >> header.version))
	header_parse_error("`Format version: " + value + "' is not a number");
      if (header.version > ConeFileHeader::current_version) {
	cerr << "This cone file is in format version " << header.version
	     << "; this LattE understands up to version "
	     << ConeFileHeader::current_version
	     << ".  It was written by a newer LattE." << endl;
	THROW_LATTE(LattException::fe_Parse, 0);
      }
      seen_version = true;
    }
    else if (key == "Ambient dimension") {
      istringstream s(value);
      if (!(s >> header.numOfVars))
	header_parse_error("`Ambient dimension: " + value + "' is not a number");
    }
    else if (key == "Homogenized") {
      if (!parse_header_bool(value, header.homogenized))
	header_parse_error("`Homogenized: " + value + "' is not 0 or 1");
    }
    else if (key == "Dualized") {
      if (!parse_header_bool(value, header.dualized))
	header_parse_error("`Dualized: " + value + "' is not 0 or 1");
    }
    else if (key == "Unbounded") {
      if (!parse_header_bool(value, header.unbounded))
	header_parse_error("`Unbounded: " + value + "' is not 0 or 1");
    }
    /* Anything else is a field of some later version: ignore it, so that
       adding a field does not by itself require a version bump. */
  }
  if (!seen_version)
    header_parse_error("no `Format version:' line");
  return header;
}

/* The input option that reads cones of the kind HEADER describes. */
static string
option_for_cones(const ConeFileHeader &header)
{
  if (!header.homogenized) {
    if (header.dualized) return "";	/* Dualized vertex cones: no option reads those. */
    return "--input-vertex-cones";
  }
  if (header.dualized) return "--input-dual-homog-cones";
  return "--input-primal-homog-cones";
}

void
checkConeFileHeader(const ConeFileHeader &header, const ConeFileHeader &expected,
		    const string &filename, const string &used_option)
{
  if (!header.present) return;	/* Written before headers existed; cannot check. */
  if (header.homogenized != expected.homogenized
      || header.dualized != expected.dualized) {
    cerr << "The cone file `" << filename << "' holds "
	 << describe_cones(header) << " (Homogenized: " << header.homogenized
	 << ", Dualized: " << header.dualized << "), but ";
    if (used_option.size()) cerr << used_option << " expects ";
    else cerr << "this input option expects ";
    cerr << describe_cones(expected) << "." << endl;
    string right_option = option_for_cones(header);
    if (right_option.size())
      cerr << "Use " << right_option << " to read this file." << endl;
    THROW_LATTE(LattException::fe_Parse, 0);
  }
  if (expected.numOfVars != 0 && header.numOfVars != 0
      && expected.numOfVars != header.numOfVars) {
    cerr << "The cone file `" << filename << "' declares ambient dimension "
	 << header.numOfVars << ", but its cones live in dimension "
	 << expected.numOfVars << "." << endl;
    THROW_LATTE(LattException::fe_Parse, 0);
  }
}

/* ----------------------------------------------------------------- */

static listVector *
readListVector(istream &in)
{
  listVector *result = NULL;
  listVector **end_p = &result;
  while (in.good()) {
    vec_ZZ v;
    skip_space(in);
    if (in.peek() != '[') break;
    in >> v;
    if (in.good()) {
      *end_p = new listVector(v);
      end_p = &(*end_p)->rest;
    }
  }
  if (result == NULL) {
    /* Nothing at all */
    return NULL;
  }
  if (result->rest == NULL
      && result->first.length() == 0) {
    /* Read [], which is meant to designate an empty list,
       rather than a list of one zero-dimensional vector. */
    freeListVector(result);
    return NULL;
  }
  return result;
}

/**
 * LattE's cone-list format, as written by printConeToFile above and read back
 * here.  An optional header block, then one block per cone, exactly:

LattE cone file.
Format version: 1
Ambient dimension: 3
Homogenized: 0
Dualized: 0
Unbounded: 0
==========
Cone.
Coefficient: 1
Vertex: [1/2 1/2 0]
Extreme rays:
[1 -1 0]
[-1 0 0]
[0 0 1]
Determinant:-1
Facets:
[]
Dual determinant:0
Lattice points in parallelepiped:
[]
==========

 * An empty list of cones is written as the single line "No cones in list."
 *
 * Vectors are bracketed and space-separated.  Ray, facet and lattice-point
 * entries are integers; only the vertex may carry fractions, written a/b, with
 * the denominator omitted when it is 1.  "[]" on its own denotes an empty
 * list, not a zero-dimensional vector.
 *
 * Both determinant fields are integers (never rationals).
 * A determinant of 0 means "not computed" or "not computable" (e.g. non-
 * simplicial cones)
 *
 * readConeFromFile consumes Coefficient, Vertex, Extreme rays, Determinant and
 * Facets.
 * Dual determinant and Lattice points in parallelepiped are written but not read
 * by readConeFromFile. 
 * subspace_generators, equalities andfacet_divisors are computed but never written.
 *
 * The header block, when present, says whether the cones are homogenized,
 * dualized or unbounded. "Ambient dimension" is Polyhedron::numOfVars, which counts the
 * homogenizing variable when the cones are homogenized.
 * The header is optional.
 *
 * Facets are outward normals with the cone on the non-positive side: the cone
 * is { x : <x,f> <= 0 for every facet f }.  They are often absent ("[]"), since
 * most producers compute rays only.
 */

listCone *
readConeFromFile(istream &in)
{
  if (!look_for(in, "Cone.")) return NULL;
  listCone *cone = createListCone();
  if (!look_for(in, "Coefficient:")) return NULL;
  in >> cone->coefficient;
  if (!in.good()) return NULL;
  if (!look_for(in, "Vertex:")) return NULL;
  /* FIXME: Actually need to handle rational data */
  skip_space(in);

  //read in a vertex.
  if (in.peek() != '[') return NULL;
  in.get(); //delete the [
  vector<RationalNTL> elements;
  RationalNTL oneTerm;

  while(in.peek() != ']')
  {
	  in >> oneTerm;
	  skip_space(in);
	  elements.push_back(oneTerm);
  }



  //vec_ZZ v;
  //in >> v;
  if (!in.good()) return NULL; //not sure what this does.--Brandon
  //ZZ denom;
  //denom = 1;
  //cone->vertex = new Vertex(new rationalVector(v, denom));

  //now save the vertex information.
  cone->vertex = new Vertex(new rationalVector(elements));

  //end of reading 1 vertex.

  if (!look_for(in, "rays:")) return NULL;
  cone->rays = readListVector(in);
  if (!look_for_capturing_determinant(in, "Facets:", cone->determinant))
    return NULL;
  cone->facets = readListVector(in);
  return cone;
}

/* ----------------------------------------------------------------- */
static void
printListConeToStream(ostream &out, listCone *cones, int numOfVars,
		      const ConeFileHeader *header)
{
  if (header) printConeFileHeaderToFile(out, *header);

  if (cones==0) out << "No cones in list.\n";

  while (cones) {
    printConeToFile(out,cones,numOfVars);
    cones = cones->rest;
  }
  out << endl;
}

void printListConeToFile(const char *fileName, listCone* cones, int numOfVars) {
  ofstream out(fileName);
  if (!out) {
    cerr << "Error opening output file `" << fileName << "' for writing in printListConeToFile!" << endl;
    exit(1);
  }
  printListConeToStream(out, cones, numOfVars, NULL);
  out.close();
  return ;
}
/* ----------------------------------------------------------------- */
void printListConeToFile(const char *fileName, listCone* cones, int numOfVars,
			 const ConeFileHeader &header) {
  ofstream out(fileName);
  if (!out) {
    cerr << "Error opening output file `" << fileName << "' for writing in printListConeToFile!" << endl;
    exit(1);
  }
  printListConeToStream(out, cones, numOfVars, &header);
  out.close();
  return ;
}
/* ----------------------------------------------------------------- */
listCone *
readListConeFromStream(istream &in, ConeFileHeader *header_out)
{
  ConeFileHeader header = readConeFileHeader(in);
  if (header_out) *header_out = header;
  listCone *result = NULL;
  listCone **tail_p = &result;
  while ((*tail_p = readConeFromFile(in)) != NULL) {
    tail_p = &(*tail_p)->rest;
  }
  return result;
}
/* ----------------------------------------------------------------- */
listCone *
readListConeFromFile(const char *filename, ConeFileHeader *header_out)
{
  ifstream in(filename);
  return readListConeFromStream(in, header_out);
}

/* ----------------------------------------------------------------- */
void
readListConeFromFile(istream &in, ConeConsumer &consumer,
		     ConeFileHeader *header_out)
{
  ConeFileHeader header = readConeFileHeader(in);
  if (header_out) *header_out = header;
  listCone *cone;
  while ((cone = readConeFromFile(in)) != NULL)
    consumer.ConsumeCone(cone);
}

/* ----------------------------------------------------------------- */
void printResidueFile(const char* fileName, listCone* cones, int numOfVars) {
  int numOfTerms;
  char outFileName[PATH_MAX];
  listVector *tmp;
  listCone *C;

  strcpy(outFileName,fileName);
  strcat(outFileName,".residue");

  ofstream out(outFileName);
  if (!out) {
    printf("Error opening output file for writing in printResidueFile!");
    exit(1);
  }
  if (cones==0) out << "No cones in list.\n";

  numOfTerms=0;

  C=cones;
  while (C) {
    numOfTerms=numOfTerms+lengthListVector(C->latticePoints);
    C=C->rest;
  }


  out << numOfVars << " " << lengthListVector(cones->rays) << " " <<
    numOfTerms << "\n\n";

  while (cones) {
    tmp=cones->latticePoints;
    while (tmp) {
      out << cones->coefficient << endl;
      printVectorToFileWithoutBrackets(out,tmp->first,numOfVars);
      printListVectorToFileWithoutBrackets(out,cones->rays,numOfVars);
      out << endl;
      tmp=tmp->rest;
    }
    cones = cones->rest;
  }
  out << endl;

  out.close();
  return ;
}

void
print_debug_vector(const vec_ZZ & v) {
   int len = v.length(); 

   cerr << "Begin vector: ["; 
   for (int i = 0; i < len; i++) {
      cerr << v[i] << ","; 
   }
   cerr << "]: End vector\n"; 
}

void
print_debug_matrix(const mat_ZZ & m) {
   int rows = m.NumRows(); 
   int cols = m.NumCols(); 

   cerr << "Begin matrix:\n"; 
   for (int i = 0; i < rows; i++) {
      cerr << "["; 
      for (int j = 0; j < cols; j++) {
         cerr << m[i][j] << ","; 
      }
      cerr << "]\n"; 
   }
   cerr << ":End matrix\n"; 
}
