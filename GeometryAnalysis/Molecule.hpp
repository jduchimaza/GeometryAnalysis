//
//  Molecule.hpp
//  GeometryAnalysis
//
//  Definition of a molecule.
//  A molecule is made up of atoms, and each atom is identified
//  by the atomic charge (zee) and has x,y,z coordinates.
//
//  Created by Juan Duchimaza on 11/1/19.
//  Copyright © 2019 chemplusplus. All rights reserved.
//

#ifndef Molecule_hpp
#define Molecule_hpp

#include <stdio.h>
#include <fstream>
using namespace std;

static double masses[] ={
    0.000000000, // void to keep index 0 null
    1.007825032, // H
    4.002603254, // He
    7.016003437, // Li
    9.012183065, // Be
    11.00930536, // B
    12.00000000, // C
    14.00307400, // N
    15.99491462, // O
    18.99840316, // F
    19.99244018  // Ne
};   // ... can fill in the rest later if needed

class Molecule
{
  public:
    int natom;
    int charge;
    double mass;
    
    int *zee;
    double **coords;
    double *com;
//    string point_group;

    void print_coords();
    double bond(int atom1, int atom2);                              // method to calculate bond length
    double angle(int atom1, int atom2, int atom3);                  // method to calculate angle
    double angle(int atom1, int atom2, int atom3, int atom4);       // method to calculate out of plane angle
    double dihedral(int atom1, int atom2, int atom3, int atom4);    // method to calculate dihedral/torsion angle
    
    
    double unitV(int cart, int atom1, int atom2);                   // helper method to give unit vector between two atoms
    void translate(double x, double y, double z);                   // helper method to translate coords by given x,y,z shift
//    void rotate(double phi);
    
 
    Molecule(string coordFile, int q);
    ~Molecule();
};
#endif /* Molecule_hpp */
