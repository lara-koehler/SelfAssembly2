//
//  SystemType.cpp
//  FrustratedSystems
//
//  Created by lara koehler on 28/09/2020.
//  Copyright © 2020 lara koehler. All rights reserved.
//

#include "SystemType.hpp"
#include "head.hpp"
#include <sstream>

#include <ctime> // For random
#include <cstdlib> // For random
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>


#include "json.hpp"
using json = nlohmann::json;

using namespace std;

void SystemType::extract_value(string &token, string const expectation){
    size_t found = (token).find(expectation);
    if (found!=std::string::npos){
        (token).erase(0, expectation.size());
        //(*token).pop_back();
    }
}

int SystemType::get_int_in_line(string token, string  const expectation){
    int result ;
    extract_value(token, expectation);
    std::stringstream stream(token);
    if (!(stream >> result)){
        cout << "Problem of conversion reading " <<expectation <<endl;
    }
    
    return result ;
    }
        


//Constructor of the class from a textfile

SystemType::SystemType() :

N_Orientations(0),
N_Edges(0),
N_Structures(0),
N_hash_structure(0),
N_particle_type(0),
Scaffold_Types(map <int, int>()),
bonds(vector<vector<int>>()),
energies(vector<float>()),
possibleFlips(0),
possibleMutations(0),
possibleStates(0)
//sequence ()
{
    }
        
void SystemType::initiate_values( string filenameStructure, string filenameLattice, string filenameEnergy, int N_particle_type_arg){
    
    std::ifstream f(filenameLattice);
    json lattice = json::parse(f);
    
    N_Edges = lattice["N_Edges"];
    N_Orientations = lattice["N_Orientations"];
    N_particle_type = N_particle_type_arg;
    N_Structures = lattice["N_Structures_"+to_string(N_particle_type)];
   
    bonds = lattice["Edges"];
    
 
      
    
// STRUCTURE MAP
fstream file(filenameStructure.c_str());
    
if (file){
    string token;
    getline(file,token);
    N_hash_structure = max(N_Orientations, N_Structures)+1;
    int counter (0);
    int counter_max(100000);
    while (token !="#END" && counter<counter_max){
        counter++;
        
        vector<string> string_list (parse(token, " "));
        if (string_list.size() !=4){
            cout << "Structure determination does not have the correct shape in file at line "<<counter<<endl;
    //                    cout<<string_list.size()<<" "<<string_list[0]<<" "<<string_list[string_list.size()-1]<<"."<<endl;
    //                    cout<<token<<endl;
        }
        else {
            int o1, o2, edge, structure ;
            o1=stoi(string_list[0]);
            o2=stoi(string_list[1]);
            edge=stoi(string_list[2]);
            structure=stoi(string_list[3]);
//            cout<<o1<<" "<<o2<<" "<<edge<<" "<<structure<<endl;
            int hash_index (hash_ints(o1, o2, edge, N_hash_structure));//o1 + o2*N_hash + edge*N_hash*N_hash);

//            cout<<o1<<" "<<o2<<" "<<edge<<" "<<structure<<" "<<hash_index<<endl;
            Scaffold_Types[hash_index]=structure;
        }
        getline(file, token);
        }
    
    if (counter>counter_max-1){
        cout <<"The structure file was not ended correctly (with #END)"<<endl;
        cout<<"Or there is too much structures, and you need to change the threshold in the code"<<endl;
    }
    }
else {
cout<<"Could not open file"<<endl;
}
        
        
//energies initialization
fstream fileEnergy(filenameEnergy.c_str());
if (fileEnergy){
    //cout <<"Energyfile found, the filename was"<<endl;
    //cout << filenameEnergy <<endl;
    
    string token;
    getline(fileEnergy,token);
    vector<string> string_list_energy (parse(token, " "));
    //cout << "string_list_energy[0]" <<string_list_energy[0]<<endl;
    int Ns = (int) string_list_energy.size()-1;
    if (Ns != N_Structures){
        cout << " PROBLEM, the lenght of the LEL vector is "<<Ns <<" and should be "<<N_Structures<<", this will cause problem in the energy computation !"<<endl;
    }
    cout<<"Energies : ";
    for (int i=0; i<Ns;i++){
        energies.push_back(stof(string_list_energy[i]));
        cout<< (energies[i])<< " ";
    }
    cout << endl;

}
else{
    cout <<"No energyfile found, the filename was"<<endl;
    cout << filenameEnergy <<endl;
}
N_Structures= (int) energies.size();
    
}

int SystemType::get_state(int particle, int orientation){
    return particle*N_Orientations +orientation; ;
}

int SystemType::get_particle(int state){
    return (int) (state-1)/N_Orientations;
}

int SystemType::get_orientation(int state){
    return (state-1)%N_Orientations +1 ;
}

void SystemType::initialize_mutation_and_flip(){
//        We build three tables :
//    possibleStates[i] is all the possible state for a particle of type i (all the orientations)
    
//    if a particle was in state i, it can be in all the states in possibleFlips[i] after a flip, and in all the states in possibleMutation[i] after a mutation.
// If it was in state 0, it means no particle, then nothing can happen
    
    possibleFlips.push_back(std::vector<int>()); //zeroth index
    possibleMutations.push_back(std::vector<int>());
    
    for (int particle=0; particle<N_particle_type; particle++){
        possibleStates.push_back(std::vector<int>());
        for (int orientation=1; orientation<N_Orientations+1; orientation++){
            possibleStates[particle].push_back(get_state(particle, orientation));
        }
    }
    
    for (int current_state=1 ; current_state<N_particle_type*N_Orientations+1;current_state++){
        
        int current_particle = get_particle(current_state) ;
        int current_orientation = get_orientation(current_state);
        
        possibleFlips.push_back(std::vector<int>()); // the index in this vector is current_state
        possibleMutations.push_back(std::vector<int>());
        
        for (int new_orientation=1; new_orientation<N_Orientations+1; new_orientation++){
            // change orientation but keep particle id
            int possible_new_state = get_state(current_particle, new_orientation);
            if (possible_new_state != current_state){
                // we add this new state unless it is the same as before
                possibleFlips[current_state].push_back(possible_new_state);
            }
        }
        
        for (int new_particle=0; new_particle<N_particle_type; new_particle++){
            // change particle but keep orientation
            int possible_new_state = get_state(new_particle, current_orientation);
            if (possible_new_state != current_state){
                // we add this new state unless it is the same as before
                possibleMutations[current_state].push_back(possible_new_state);
            }
        }
    }
    
//    cout<<"Verify flips and mutation"<<endl;
//    for (int myState=1 ; myState<possibleMutations.size() ; myState++){
//        int myParticle=get_particle(myState);
//        int myOrientation=get_orientation(myState);
//        cout<<"State "<<myState<<"("<<myParticle<<" "<<myOrientation<<")"<<endl;
//        cout<<"Possible flips : ";
//        for (int iflip=0; iflip<possibleFlips[myState].size(); iflip++){
//            cout<<possibleFlips[myState][iflip]<<" ";
//        }
//        cout<<endl;
//        cout<<"Possible mutations : ";
//        for (int imut=0; imut<possibleMutations[myState].size(); imut++){
//            cout<<possibleMutations[myState][imut]<<" ";
//        }
//        cout<<endl;
//    }

    
}

int SystemType::possible_new_state(int current_state, bool isMutation){
    // if isMutation, we do a mutation : particle type changes but orientation stays the same
    // if not isMutation, we do a flip : particle type stays the same but orientation changes
    
    //Range of the random number we will draw
    int N_possibilities ;
    if (isMutation){
        N_possibilities = N_particle_type-1;
    }
    else {
        N_possibilities = N_Orientations-1;
    }
    
    // Choose the new state in the tables possibleFlips or possibleMutations
    int new_state ;
    int id_new_state = rand()%N_possibilities;
    if (isMutation){
        new_state = possibleMutations[current_state][id_new_state];
    }
    else{
        new_state = possibleFlips[current_state][id_new_state];
    }
    
    return new_state ;
}

int SystemType::random_state(int particle){
    // knowing the type of particle, choose a random orientation and returns the correponding site state
    int N_possibilities = N_Orientations;
    int id_orientation = rand()%N_possibilities;
    return possibleStates[particle][id_orientation];
}



float SystemType::get_energy(int index) const{
    if  (index>=N_Structures){
        cout <<"problem you ask for energy of index "<<index<<endl;
    }
    return energies[index];
}

vector<float> SystemType::get_energies() const {
    return energies ;
}

int SystemType::get_scaffold_value(int const o1, int const o2, int const edge) const{
    int hash_index =hash_ints(o1, o2, edge, N_hash_structure);
    int test = (int) Scaffold_Types.count(hash_index);
    if (test ==0){
        cout << "wrong requested scaffold"<<endl;
        cout <<"o1="<<o1<<" o2="<<o2<<" edge="<<edge<<endl;
        cout<< "hash index was "<<hash_index;
    }
//    cout<<"get_scaffold "<<o1<<" "<<o2<<" "<<edge<<" "<<Scaffold_Types.at(hash_index)<<endl;
    return (Scaffold_Types.at(hash_index));
    
    /*if (result ==0 &&	 (o1!=0 || o2!=0 || edge!=0)){
        cout << "The combination of orientation and scaffold asked is not possible for "<<o1<<" "<<o2<<" "<<edge<<endl;
    }
    return result;*/
}



