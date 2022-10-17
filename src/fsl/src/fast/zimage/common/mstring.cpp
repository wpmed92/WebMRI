/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include <fstream>
#include <strstream>
#include <cstdio>
#include "mstring.h"
using namespace std;

const char* FileName(const char* fname)
{
	const char* rname = strrchr(fname, '/');
	if(rname == 0) rname = strrchr(fname, '\\');
	if(rname == 0) rname = fname;
	else rname++;
	return rname;
}

const char* ExtName(const char* fname)
{
	const char* rname = strrchr(fname, '.');
	if(rname == 0) rname = fname+strlen(fname);
	else rname++;
	return rname;
}

string BaseName(const char* fname)
{
	string name = fname;
	const char* rname = strrchr(fname, '.');
	if(rname == 0) return name;
	name.assign(fname, rname);
	return name;
}

string Path(const char* fname)
{
	string path;
	const char* rname = strrchr(fname, '/');
	if(rname == 0) rname = strrchr(fname, '\\');
	if(rname != 0) path.assign(fname, rname+1);
	return path;
}

bool GetParameterString(const char* filename, const char* keyname, char* string)
{
	fstream		file(filename, ios::in);
	char		line[200], *con;
	int			keylength = strlen(keyname);

	if(file.fail()) return false;

	for(file.getline(line,200); !file.eof(); file.getline(line,200))
	{
		if(line[0] == '#') continue;

		if(!strncmp(keyname, line, keylength)) break;
	}

	bool ret = false;
	if(!file.eof())
	{
		if((con = strchr(line, '=')) != NULL || (con = strchr(line, ':'))!=NULL)
		{
			con++;

			for(;*con == ' ' || *con == '\t'; con++);
			if(con-line < (int)strlen(line))
			{
				ret = true;
				strcpy(string, con);
			}
		}
	}

	file.close();

	return ret;
}

bool GetParameterNumber(const char* filename, const char* keyname, int& number)
{
	fstream		file(filename, ios::in);
	char		line[200], *con;
	int			keylength = strlen(keyname);

	if(file.fail()) return false;

	for(file.getline(line,200); !file.eof(); file.getline(line,200))
	{
		if(line[0] == '#' || !line[0]) continue;

		if(!strncmp(keyname, line, keylength)) break;
	}

	bool ret = false;
	if(!file.eof())
	{
		if((con = strchr(line, '=')) != NULL || (con = strchr(line, ':'))!=NULL)
		{
			ret = true;
			sscanf(con+1, "%d", &number);
		}
	}

	file.close();

	return ret;
}

bool GetParameterNumber(const char* filename, const char* keyname, float& number)
{
	fstream		file(filename, ios::in);
	char		line[200], *con;
	int			keylength = strlen(keyname);

	if(file.fail()) return false;

	for(file.getline(line,200); !file.eof(); file.getline(line,200))
	{
		if(line[0] == '#') continue;

		if(!strncmp(keyname, line, keylength)) break;
	}

	bool ret = false;
	if(!file.eof())
	{
		if((con = strchr(line, '=')) != NULL || (con = strchr(line, ':'))!=NULL)
		{
			ret = true;
			sscanf(con+1, "%f", &number);
		}
	}

	file.close();

	return ret;
}

bool GetParameterNumber(const char* filename, const char* keyname, double& number)
{
	fstream		file(filename, ios::in);
	char		line[200], *con;
	int			keylength = strlen(keyname);

	if(file.fail()) return false;

	for(file.getline(line,200); !file.eof(); file.getline(line,200))
	{
		if(line[0] == '#') continue;

		if(!strncmp(keyname, line, keylength)) break;
	}

	bool ret = false;
	if(!file.eof())
	{
		if((con = strchr(line, '=')) != NULL || (con = strchr(line, ':'))!=NULL)
		{
			ret = true;
			sscanf(con+1, "%lf", &number);
		}
	}

	file.close();

	return ret;
}

bool isNumber(const string& keyword) 
{
	strstream tmp; tmp << keyword;
	float number;
	if(!(tmp >> number)) return false; 
	return true;
}

bool isIdentifier(const string& keyword) 
{
	if (!isalpha(keyword[0])) return false;
	for (unsigned int count = 1; count < keyword.length(); count++) 
	{
		if ((!isalnum(keyword[count])) && (keyword[count]!='_')) return false;
	}
	
	return true;
}

string NextToken(fstream& f_in)
{
	char buf[1024];
	
	char c; string token;
	int lastchar = 1;
	if (!f_in || f_in.eof()) return(token); 
	while (f_in.get(c)) 
	{
		if (f_in.eof()) return(token);
		switch (c) 
		{
			case '\n': 
			case ' ' :
			case '\t':
			case '\a':
			case '\b':
			case '\v':
			case '\f':
			case ',':   
			case ';': 
			case '\r': if ( lastchar == 0 ) return(token); break;
			case '#': f_in.getline(buf, 1024); break;
			case '/': if(f_in.peek() == '/') f_in.getline(buf, 1024); break;
			case '[':   
			case ']':   
			case '{':   
			case '}':  token = c; return(token);
			case '!':  while ( c != '\n' && f_in.get(c)); break;
			default: 
			{
				token = token + c;	lastchar = 0;
				char cc =  f_in.peek();
				if(!isalnum(cc) && cc!='_')
				{
					if ( lastchar == 0 ) return(token);
					else { f_in.get(c); token = c;  return(token); }
				}
			}
		}
	}
	
	return(" ");
}
