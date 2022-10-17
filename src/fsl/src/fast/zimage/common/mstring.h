/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file declares string functions.

#ifndef __MSTRING_H__
#define __MSTRING_H__

#include <string>
#include <fstream>

bool isNumber(const std::string& keyword);
bool isIdentifier(const std::string& keyword);
std::string NextToken(std::fstream& f_in);

const char* FileName(const char* fname);
const char* ExtName(const char* fname);
std::string BaseName(const char* fname);
std::string Path(const char* fname);

bool GetParameterString(const char* filename, const char* keyname, char* string);
bool GetParameterNumber(const char* filename, const char* keyname, int& number);
bool GetParameterNumber(const char* filename, const char* keyname, float& number);
bool GetParameterNumber(const char* filename, const char* keyname, double& number);

#endif
