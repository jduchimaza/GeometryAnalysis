//
//  Molecule.cpp
//  GeometryAnalysis
//
//  Created by Juan Duchimaza on 11/1/19.
//  Copyright © 2019 chemplusplus. All rights reserved.
//

#include "Molecule.hpp"
#include <cstdio>
#include <fstream>
#include <math.h>

void Molecule::print_coords(){
    printf("Z\tx\t\t y \t\t z\n");
    for(int i=0; i < natom; i++){
      printf("%d %20.12f %20.12f %20.12f\n", zee[i], coords[i][0], coords[i][1], coords[i][2]);
    }
}

double Molecule::bond(int atom1, int atom2){
    /*
     Bond length = sqrt(x^2 + y^2 + z^2)
     where x = x1 - x2 (same for y, same for z)
     */
    double x2 = pow(coords[atom1][0] - coords[atom2][0],2);
    double y2 = pow(coords[atom1][1] - coords[atom2][1],2);
    double z2 = pow(coords[atom1][2] - coords[atom2][2],2);
    return sqrt(x2+y2+z2);
}

double Molecule::angle(int atom1, int atom2, int atom3){
    /*
     Bond angle: have three atoms i-j-k, where j is central atom
     angle is defined as cos(theta,ijk)=vec(j->i)•vec(j->k)
     */
    double dotprod = unitV(0, atom2, atom1)*unitV(0, atom2, atom3) + unitV(1, atom2, atom1)*unitV(1, atom2, atom3) + unitV(2, atom2, atom1)*unitV(2, atom2, atom3);
    double angl = acos(dotprod); // radians
    angl = angl*(180.0/acos(-1.0)); // degrees
    return angl;
}

double Molecule::angle(int atom1, int atom2, int atom3, int atom4){
    // Overloaded angle method - should give "out of plane angle", i.e. the angle of an atom i from a plane formed by atoms jkl
    // sin(theta,ijkl) = [(vec(k->j) X vec(k->l))/sin(theta,jkl)] • vec(k->i)
    // so need sin(angle(j,k,l))
    // unitV for the angles
    // maybe I should define a dot-prod and cross-prod function
    // use asin on the whole RHS of the expression
    //
    double xprod[3]; //the cross product vector
    xprod[0] = (unitV(1, atom3, atom2)*unitV(2,atom3,atom4))-(unitV(2,atom3,atom2)*unitV(1,atom3,atom4));
    xprod[1] = (unitV(2, atom3, atom2)*unitV(0,atom3,atom4))-(unitV(0,atom3,atom2)*unitV(2,atom3,atom4));
    xprod[2] = (unitV(0, atom3, atom2)*unitV(1,atom3,atom4))-(unitV(1,atom3,atom2)*unitV(0,atom3,atom4));
    double dotprod=xprod[0]*unitV(0,atom3,atom1) + xprod[1]*unitV(1,atom3,atom1) + xprod[2]*unitV(2,atom3,atom1);
    double angl = dotprod/sin(angle(atom2,atom3,atom4)/(180.0/acos(-1.0)));
    
    if(angl < -1.0) angl = asin(-1.0);
    else if(angl > 1.0) angl = asin(1.0);
    else angl = asin(angl);
    
    angl = angl*(180.0/acos(-1.0)); // turn radians to degrees
    return angl;
}

double Molecule::unitV(int cart, int atom1, int atom2){
    return -(coords[atom1][cart] - coords[atom2][cart])/bond(atom1,atom2);
}

// Constructor
Molecule::Molecule(string coordFile, int q){
    charge = q;
    // Open the file
    std::ifstream inFile(coordFile);
    assert(inFile.good());
    
    // Read in the number of atoms
    inFile >> natom;
    
    // Allocate space for vars
    zee = new int[natom];
    coords = new double* [natom];
    
    for(int i=0; i<natom; i++){
        coords[i] = new double[3];
        // Read in: zee     x             y             z
        inFile >> zee[i] >> coords[i][0] >> coords[i][1] >> coords[i][2];
    }
    
    // Close the file
    inFile.close();
}

// Destructor
Molecule::~Molecule()
{
  delete[] zee;
  for(int i=0; i < natom; i++)
    delete[] coords[i];
  delete[] coords;
}
