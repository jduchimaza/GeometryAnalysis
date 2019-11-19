//
//  main.cpp
//  GeometryAnalysis
//
//  Created by Juan Duchimaza on 11/1/19.
//  Copyright © 2019 chemplusplus. All rights reserved.
//
//  This is an implementation of the Programming Project from
//  http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project1
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include "Molecule.hpp"
#include <math.h>

using namespace std;
int main(int argc, const char * argv[]) {
    // Print some header
    cout << "Geometry Analysis\n";

    // Define my working directory
    string pwd = "/Users/juanduchimaza/Documents/git/coding/GeometryAnalysis/GeometryAnalysis/";
    // Define the file I'll be using (maybe this can be passed to argv)
    string inFile = pwd+"acetaldehyde.xyz";
    
    // Initialize acetaldehyde molecule instance
    Molecule acet(inFile, 0);
    
    // Print the coordinates
    acet.print_coords();
     
    // Print bond lengths
    cout << "\nInteratomic distances: \n";
    for(int i = 0; i < acet.natom; i++){
        for(int j = 0; j < i; j++){
            printf("Atoms %2d and %2d: %8.5f\n", i, j, acet.bond(i,j));
            //cout << "Interatomic distance between atoms " << i << " and " << j << ": " << acet.bond(i,j) << "\n";
        }
    }
    
    // Print bond angles
    cout << "\nBond Angles:\n";
    for(int i = 0; i < acet.natom; i++){
        for(int j = 0; j < i; j++){
            for(int k=0; k < j; k++){
                if(acet.bond(i,j) < 4.0 && acet.bond(j,k) < 4.0)
                    printf("Atoms %2d, %2d, %2d = %10.6f\n",i,j,k,acet.angle(i, j, k)*(180.0/acos(-1.0)));
            }
        }
    }
    
    // Print out of plane angles
    cout << "\nOut of plane Angles:\n";
    for(int i = 0; i < acet.natom; i++){
        for(int j = 0; j < acet.natom; j++){
            for(int k = 0; k < acet.natom; k++){
                for(int l = 0; l < j; l++){
                    if(i!=j && i!=k && i!=l && j!=k && k!=l && acet.bond(i,k) < 4.0 && acet.bond(k,j) < 4.0 && acet.bond(k,l) < 4.0)
                        printf("Atoms %2d, %2d, %2d, %2d = %10.6f\n",i,j,k,l,acet.angle(i, j, k, l)*(180.0/acos(-1.0)));
                }
            }
        }
    }
    
    // Print dihedral/torsion angles
    cout << "\nDihedrals:\n";
    for(int i=0; i < acet.natom; i++){
        for(int j=0; j < i; j++){
            for(int k=0; k < j; k++){
                for(int l=0; l < k; l++){
                    if(acet.bond(i,j) < 4.0 && acet.bond(j,k) < 4.0 && acet.bond(k,l) < 4.0){
                        printf("Atoms %2d, %2d, %2d, %2d = %10.6f\n",i,j,k,l,acet.dihedral(i, j, k, l)*(180.0/acos(-1.0)));
                    }
                }
            }
        }
    }

    cout << "\nMolecular Mass:\n";
    printf("%10.10f\n",acet.mass);
    
    cout << "\nCenter of Mass:\n";
    printf("x: %12.8f, y: %12.8f, z: %12.8f\n", acet.com[0], acet.com[1], acet.com[2]);
    acet.translate(acet.com[0], acet.com[1], acet.com[2]);
    
    cout << "\nNew Translated Coords:\n";
    acet.print_coords();
    
    printf("\n\n");
    return 0;
}
