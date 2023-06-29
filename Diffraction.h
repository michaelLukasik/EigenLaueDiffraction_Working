#ifndef DIFFRACTION_H
#define DIFFRACTION_H

#include <string>
#include <vector>
#include <Eigen/Dense>


bool debug;
std::string cellName;
std::string cellType;
std::string cellCentering;

// Unit cells in x, y, z directions 
int nx;
int ny;
int nz;

double d;// Interatomic Distance (typically in angstroms i.e. E-10)
double wallXPosition;// wall at x = wallXPos
double dzdy;// dzdy->step size for wall spacing
double wallLength; // wall is wallLen x wallLen meters
double lambda; // wavelength (meters)
double A;  // Amplitude (arbitrary units)
double k; // wavenumber (1/m)
double omega; //
double theta;
double phi;
double psi;

int wallDivisions;

// Matrix dealings

Eigen::MatrixXd screen;

#endif