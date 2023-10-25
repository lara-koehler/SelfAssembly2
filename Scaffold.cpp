//
//  Scaffold.cpp
//  FrustratedSystems
//
//  Created by lara koehler on 23/09/2020.
//  Copyright © 2020 lara koehler. All rights reserved.
//

#include "Scaffold.hpp"
#include "head.hpp"
#include "SystemType.hpp"
#include <vector>
#include <algorithm>
using namespace std ;




Scaffold::Scaffold(int edgeRef, SystemType typeOfSystem) : m_edgeRef(edgeRef), m_structureRef(0), m_orientationList(vector <int> (2,0)), typeOfSystem(typeOfSystem){
    if (m_edgeRef!=0){
        m_structureRef=typeOfSystem.get_scaffold_value(0, 0, m_edgeRef);
    }
        
}


int Scaffold:: getStructure() const
{
    return m_structureRef;
}

int Scaffold::recompute_and_getStructure() {
    m_structureRef = typeOfSystem.get_scaffold_value(m_orientationList[0], m_orientationList[1], m_edgeRef);
    return m_structureRef;
}


void Scaffold::actualize_scaffold(int newOrientation, int newIndex){
    this->change_orientations(newOrientation, newIndex);
    this->actualize_structure();
    
}


void Scaffold::change_orientations(int newOrientation, int newIndex){
    if (newIndex !=0 && newIndex !=1){
        cout << "The index given to actualize scaffold was not correct"<<endl;
    }
    if (m_edgeRef==0){
    m_orientationList[0]=newOrientation;
    m_orientationList[1]=newOrientation;
    }
    else{
    m_orientationList[newIndex]=newOrientation;
    }
    
}

void Scaffold::actualize_structure(){
    m_structureRef = typeOfSystem.get_scaffold_value(m_orientationList[0], m_orientationList[1], m_edgeRef);
}



float Scaffold::get_energy() const {
//    if (m_edgeRef>0){
//        bool bool1 = m_orientationList[0]>0;
//        bool bool2 = m_orientationList[1]>0;
//        cout << "or " << bool1 <<" "<<bool2	<<" st ";
//        cout <<m_structureRef <<" ";
//        cout<<" energy "<< typeOfSystem.get_energy(m_structureRef)<<endl;
//        
//    }
    return typeOfSystem.get_energy(m_structureRef);
}



void Scaffold::printAttributes() const {
    cout << "Adress : "<<this;
    cout << ", EdgeType : "<<m_edgeRef;
    vector <int> :: const_iterator it;
    cout <<", Orientations of member sites : ";
    for (it = m_orientationList.begin() ; it!=m_orientationList.end() ;  it++){
        cout  <<(*it)<< " ";
    }
    cout << ", Resulting Structure : "<<m_structureRef;
    cout <<endl;
}



