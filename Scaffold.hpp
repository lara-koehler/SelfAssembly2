//
//  Scaffold.hpp
//  FrustratedSystems
//
//  Created by lara koehler on 23/09/2020.
//  Copyright © 2020 lara koehler. All rights reserved.
//

#ifndef Scaffold_hpp
#define Scaffold_hpp

#include <stdio.h>
#include <string>
#include <vector>
//#include "Site.hpp"
#include "SystemType.hpp"


class Scaffold {
    
    public :
    //Scaffold();
    
    Scaffold(int edgeRef, SystemType typeOfSystem) ;
    //Scaffold(int m_edgeRef , int structureRef, std::vector <int> orientationList);
    
    //Scaffold(/*std::vector <Site*> listOfSites*/);
    
    int getStructure() const ;
    int recompute_and_getStructure() ;
    
    float get_energy() const ;
    
    void printAttributes() const ;
    
    //static void printStructureRefs(std::vector <Scaffold *> listToPrint);
    
    
    //void initialize_1p(int orientation);
    //void initialize_2p(int orientation1, int index1, int orientation2, int index2);
    
    
    private :
    int const m_edgeRef;
    int m_structureRef; //should be removed and put energy instead. Still, we can have a function that access the 
    std::vector <int> m_orientationList ;
    //std::vector <Site*> m_memberSites; //points to the sites that are member of the scaffold
    const SystemType typeOfSystem ;
    
    
    
    void actualize_scaffold(int newOrientation, int newIndex);//changes the value of orientation of one of its member site, and recomputes its new structure.
    void actualize_structure();
    void change_orientations(int newOrientation, int newIndex);
    
    
    friend class Site;
    
};



#endif /* Scaffold_hpp */
