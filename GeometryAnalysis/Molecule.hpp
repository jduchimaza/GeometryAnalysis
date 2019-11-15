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

class Molecule
{
  public:
    int natom;
    int charge;
    int *zee;
    double **coords;
//    string point_group;

    void print_coords();
    double bond(int atom1, int atom2);
    double angle(int atom1, int atom2, int atom3);
    double angle(int atom1, int atom2, int atom3, int atom4);
    
    double unitV(int cart, int atom1, int atom2);
    
    
//    void rotate(double phi);
//    void translate(double x, double y, double z);
//    double torsion(int atom1, int atom2, int atom3, int atom4);
 
    Molecule(string coordFile, int q);
    ~Molecule();
};
#endif /* Molecule_hpp */
