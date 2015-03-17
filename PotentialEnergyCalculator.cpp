// Reads in an MDR1 halo data file in order to calculate its potential energy.
// August 20, 2013
// Austen Groener

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cstring>
#include <cmath>
#include <cstdlib>
using namespace std;


vector<string> split(const string &s,char delim) 
{
  stringstream ss(s);
  string item;
  vector<string> elems;
  while(getline(ss,item,delim))
    {
      elems.push_back(item);
    }
  return elems;
}


int main(int argc,char* argv[])
{
  char cword[25];
  ifstream textfile;
  vector< vector<string> > data;
  char* skiplines = "particleId,snapnum,x,y,z,vx,vy,vz,phkey";
  textfile.open(argv[1]);
  char tmp[1024];
  while(!textfile.eof())
    {
      stringstream line;
      textfile.getline(tmp,1024);

      if (tmp[0] == '#')
	{
	  continue;
	}

      if (!strcmp(tmp,skiplines))
	{
	  continue;
	}

      line << tmp;
      vector<string> row=split(line.str(),',');
      data.push_back(row);
    }

  int colstart = 2;
  int colend = 4;
  double PE = 0;
  double p_mass = 8.721e9 * 1.989e30; // mass of sim particle in kg
  double mpc_to_m = 3.08567758e22; // meters per mpc
  double G = 6.67384e-11; // grav. constant in mks units
  double prefactor = (-1 * G * pow(p_mass,2)); // prefactor for quicker multiplication
  
  // Calculating PE
  for(int i = 0; i < data.size()-1; i++) // -1 is for empty line at the end of every file
    {

      double x_i = (double)atof(data[i][2].c_str());
      double y_i = (double)atof(data[i][3].c_str());
      double z_i = (double)atof(data[i][4].c_str());
       
      for(int j = 0; j < i; j++)
	{

	  double x_j = (double)atof(data[j][2].c_str());
	  double y_j = (double)atof(data[j][3].c_str());
	  double z_j = (double)atof(data[j][4].c_str());

	  double dist_betw = sqrt((x_i-x_j)*(x_i-x_j) + (y_i-y_j)*(y_i-y_j) + (z_i-z_j)*(z_i-z_j)) * mpc_to_m;
	  
	  PE += prefactor / (dist_betw);  
  
	}
    }
  
  cout << PE << " " << argv[1] ;
  textfile.close();
  return 0;
}


