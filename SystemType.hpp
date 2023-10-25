//
//  SystemType.hpp
//  FrustratedSystems
//
//  Created by lara koehler on 28/09/2020.
//  Copyright © 2020 lara koehler. All rights reserved.
//

#ifndef SystemType_hpp
#define SystemType_hpp
#include <vector>
#include <stdio.h>
#include <map>
#include <string>
#include "head.hpp"



class SystemType {

    private :
    int N_Orientations;
    int N_Edges;
    int N_Structures;
    int N_hash_structure;
    int N_particle_type ;
    std::map <int, int>  Scaffold_Types;
    std::vector<std::vector <int> > bonds ;
    std::vector<float> energies;
//    const std::string sequence ;
    
    //static const std::string filename ;
    //static const std::string filenamelattice ;

    public :
    SystemType();
    
    void initiate_values( std::string filenameStructure, std::string filenameLattice, std::string filenameEnergy, int N_particle_type_arg);
    
    float get_energy(int index) const;
    
    std::vector<float> get_energies() const;
    
    int get_N_Structures() const{
        return N_Structures ;
    };
    
    int get_state(int particle, int orientation);
    int get_particle(int state);
    int get_orientation(int state);
    
    int possible_new_state(int current_state, bool isMutation) ;
    int random_state(int particle) ;
    void initialize_mutation_and_flip();
    
//    const std::string get_sequence() const{
//        return sequence ;
//    }
    
    std::vector<std::vector<int> > get_bonds() const {
        return bonds;
    };
    
    int get_scaffold_value(int const o1, int const o2, int const edge) const;

    //int get_2particles_scaffold (int orientation1, int orientation2, int edge);
        
    int getNOrientation() const {
        return N_Orientations;
    }
    
    static void extract_value(std::string &token, std::string expectation);
    // given a line of characters, return what is behind the expectation string
    
    static int get_int_in_line(std::string token, std::string  const expectation);

    //static std::vector<std::string> parse(std::string s, std::string delimiter);
    
    
    std::vector<std::vector<int>> possibleFlips ; //for a given state i, possibleFlips[i] is all the other state particle can have after flip
    std::vector<std::vector<int>> possibleMutations ; //for a given state i, possibleFlips[i] is all the other state particle can have after mutation
    
    std::vector<std::vector<int>> possibleStates; // all the possible states given the particle id
};



#endif /* SystemType_hpp */
