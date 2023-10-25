//
//  Site.cpp
//  FrustratedSystems
//
//  Created by lara koehler on 22/09/2020.
//  Copyright © 2020 lara koehler. All rights reserved.
//

#include "Site.hpp"
#include "Scaffold.hpp"
#include "head.hpp"
#include <ctime> // For random
#include <cstdlib> // For random
#include <iostream>
#include <algorithm>
#include <cmath>
#include <time.h>

//#include "SystemType.hpp"
using namespace std ;
class Scaffold ;

/*Site::Site(int orientationRef, vector<vector <Site*>> neighbourList) : m_orientationRef (orientationRef), m_neighbourList (vector<vector <Site *>>(0)), hashedPosition(0), belongingScaffolds(vector <Scaffold*> (0)), indexInBelongingScaffold(vector <int>()), belongingExtendedScaffolds(vector <Scaffold*> (0)), indexInBelongingExtendedScaffold(vector <int>())
{
    vector <Site *> :: iterator i;
    for (i=neighbourList.begin(); i != neighbourList.end() ; ++i){
        (*this).setNeighbourRelation(**i);
    }
}*/

Site::Site() : m_orientationRef(0), m_neighbourList(vector<vector <Site *>> (0)), hashedPosition(0), belongingScaffolds(vector<vector <Scaffold*>> (0)),  indexInBelongingScaffold(vector<vector <int>>())
{
    for (int i=0;i<maxDegreeOfNeighbourhood;i++){
        m_neighbourList.push_back(vector <Site*>(0));
        belongingScaffolds.push_back(vector<Scaffold*>(0));
        indexInBelongingScaffold.push_back(vector <int>());
    }
//    cout<< this << " has "<< m_neighbourList.size()<<" type of neighbors"<<endl;

}


vector <Site*> Site::getNeighbours (int degree){
    return m_neighbourList[degree];
}


/*
Site*  Site::getAdress() const
{
    return this;
}
 */

void Site::setScaffoldAppartenance(Scaffold* s_adress, int indexInS, int degreeOfNeighbourhood){
    belongingScaffolds[degreeOfNeighbourhood].push_back(s_adress);
    indexInBelongingScaffold[degreeOfNeighbourhood].push_back(indexInS);
} 

  
                                  
                    
void Site::printAttributes(int const Nhash) const{
    cout << endl<< "Site description"<<endl;
    cout << "The adress of the site is : " <<this <<endl;
    cout << "The orienation is : " <<m_orientationRef <<endl;
    
    cout << "The hashed position is : "<<hashedPosition<<endl;
    vector <int> real_pos (this->get_real_space_position(Nhash));
    cout << "The real position is : "<<real_pos[0]<<" "<<real_pos[1]<<" "<<real_pos[2]<<endl;
    
    //Site::printVectorOfSite(m_neighbourList);
    cout << "Neighbours and their orientation : "<<endl;
    vector <Site *> ::const_iterator iter;
    for (iter = m_neighbourList[0].begin() ; iter!=m_neighbourList[0].end() ; iter++){
        cout << *iter << " -> "<< (**iter).m_orientationRef <<endl;
    }
    
    cout << "Second neighbours and their orientation : "<<endl;
    vector <Site *> ::const_iterator iter2;
    for (iter2 = m_neighbourList[1].begin() ; iter2!=m_neighbourList[0].end() ; iter2++){
        cout << *iter2 << " -> "<< (**iter2).m_orientationRef <<endl;
    }
    //printVectorOfAdress(m_neighbourList);
    
    cout << "Scaffolds and their structures : "<<endl;
    vector <Scaffold *> ::const_iterator it;
    for (it = belongingScaffolds[0].begin() ; it!=belongingScaffolds[0].end() ; it++){
        cout << *it << " -> "<< (*it)->getStructure() <<endl;
    }
    //printVectorOfAdress(belongingScaffolds);
    //Scaffold::printStructureRefs(belongingScaffolds);
    
    /*vector <Site *> :: const_iterator i ;
    for (i = m_neighbourList.begin(); i != m_neighbourList.end() ; ++i){
        cout << *i <<endl;}
    cout <<endl;*/
}

/*
void Site::printVectorOfSite(vector <Site *> listToPrint){
    //if listToPrint.site()!=0){
    vector <Site *> :: const_iterator i ;
    for (i = listToPrint.begin(); i != listToPrint.end() ; ++i){
        cout << *i <<endl;}
    cout <<endl;
        
}*/


bool Site::setNeighbourRelation(Site& neighbourSite, int degreeOfNeighbour){
    bool madeBond = false;
    bool madeReverseBond = false;
    if (!(this==&neighbourSite)){
    //vector <Site*> &listNeighbour (neighbourSite.m_neighbourList);
    if (! ( (*this).isNeighbourOf(neighbourSite,degreeOfNeighbour) ))
       { //if the current site is not already a neighbour of neighbourSite
            neighbourSite.m_neighbourList[degreeOfNeighbour].push_back(this); //add this to the list of neighbour adress for the neighbourSite
            madeBond = true;
    }
    if (! (neighbourSite.isNeighbourOf(*this, degreeOfNeighbour)) )
        { //if the neighbourSite is not already a neighbour of the current site
        m_neighbourList[degreeOfNeighbour].push_back(&neighbourSite); //add the neighbour adresss to the list of the current site neighbours adresses
            madeReverseBond=true;
    }
    }
    //else{cout << "found equality while setting neighbour relation"<<endl;}
    return(madeBond && madeReverseBond);
}




void Site::set_position(int const i, int const  j,int  const  k, int const  N){
    hashedPosition=hash_ints(i, j, k, N);
    
}

vector<int> Site::get_real_space_position(int const N) const {
    return(unhash_ints(hashedPosition, N));
}

//void Site::fillWithRandomValue(int numberParticleType){
//    setFullSite( ( rand()%(numberParticleType) ) +1 );
//}
//
//void Site::fillWithNewRandomValue (int numberParticleType, int currentOrientation){
//    setFullSite(  1 + ( currentOrientation + rand()%(numberParticleType-1) )%numberParticleType );
//}

                                                                                                                       



void Site::setFullSite(int orientationRef){
//    int oldOrientation = m_orientationRef;
   // cout << m_orientationRef <<"->"<<orientationRef<<endl;;
    m_orientationRef = orientationRef ;

//    this->newFullSite();
    for (int degreeOfNeighbourhood=0 ; degreeOfNeighbourhood<maxDegreeOfNeighbourhood ; degreeOfNeighbourhood++){

        vector <Scaffold *> :: iterator it;
        vector <int>:: iterator indexIterator;
        indexIterator=indexInBelongingScaffold[degreeOfNeighbourhood].begin();
        for (it = belongingScaffolds[degreeOfNeighbourhood].begin(); it!= belongingScaffolds[degreeOfNeighbourhood].end(); it++){
//            int oldStructure =(*it)->m_structureRef;
            
           // cout << (*it)->m_structureRef<<"->";
            //(*it)->actualize_scaffold(m_orientationRef,indexInBelongingScaffold[degreeOfNeighbourhood].at(*indexIterator));
            (*it)->actualize_scaffold(m_orientationRef,*indexIterator);
//            int newStructure=(*it)->m_structureRef;
            indexIterator++;
        }
        
}
}

void Site::changeScaffoldOrientation(){

     for (int degreeOfNeighbourhood=0 ; degreeOfNeighbourhood<maxDegreeOfNeighbourhood ; degreeOfNeighbourhood++){

        vector <Scaffold *> :: iterator it;
//        vector <int>:: iterator indexIterator;
//        indexIterator=indexInBelongingScaffold[degreeOfNeighbourhood].begin();
        int indexIterator = 0;
        for (it = belongingScaffolds[degreeOfNeighbourhood].begin(); it!= belongingScaffolds[degreeOfNeighbourhood].end(); it++){
            
//            (*it)->change_orientations(m_orientationRef,indexInBelongingScaffold[degreeOfNeighbourhood].at(*indexIterator));
            (*it)->change_orientations(m_orientationRef,indexInBelongingScaffold[degreeOfNeighbourhood][indexIterator]);
            indexIterator++;
        
        }
}
}

void Site::changeScaffoldStructure(){
//    if (m_orientationRef==0){
//          cout<<"New empty site ";
//      }
     for (int degreeOfNeighbourhood=0 ; degreeOfNeighbourhood<maxDegreeOfNeighbourhood ; degreeOfNeighbourhood++){

        vector <Scaffold *> :: iterator it;
        for (it = belongingScaffolds[degreeOfNeighbourhood].begin(); it!= belongingScaffolds[degreeOfNeighbourhood].end(); it++){
//            if (m_orientationRef==0){
//             cout << (*it)->getStructure()<<"->";
//            }
            
            (*it)->actualize_structure();

//            if (m_orientationRef==0){
//             cout << (*it)->getStructure()<<" ";
//            }
        }
}
}





int Site::getOrientationRef() const {
    return m_orientationRef;
}

int Site::getSiteStateTemporary() {
    return m_orientationRef;
}

float Site::getSiteEnergy(){
    float energy =0;
    vector <Scaffold *> :: iterator it;
    
//    cout << "Site energy ";
    for (it = belongingScaffolds[0].begin(); it!= belongingScaffolds[0].end(); it++){
        energy+=(*it)->get_energy();
//        cout<<energy<<", ";
    
    }
//    cout<<endl;
    //cout<<"Etot "<<energy<<endl;
    return energy;
}

void Site::setOrientation(int orientationValue){
    m_orientationRef=orientationValue;
}

void Site::putParticleIn(Site &destinationSite, int new_state){
    // the destination site has the new_state
    // the original site has the state of the destination site (usually, empty)

    
    int orientationCopy = destinationSite.m_orientationRef;
    //Change both orientation at the same time
    //destinationSite.setOrientation(m_orientationRef);
    destinationSite.setOrientation(new_state);
    
    m_orientationRef = orientationCopy ;
    
    destinationSite.changeScaffoldOrientation();
    this->changeScaffoldOrientation();
    
    destinationSite.changeScaffoldStructure();
    this->changeScaffoldStructure();
    
//    for (int degreeOfNeighbourhood=0 ; degreeOfNeighbourhood<maxDegreeOfNeighbourhood ; degreeOfNeighbourhood++){
//            vector <Scaffold *> :: iterator it;
//
//            //Actualize the scaffolds of Site
//            vector <int>:: iterator indexIterator;
//            indexIterator=indexInBelongingScaffold[degreeOfNeighbourhood].begin();
//            for (it = belongingScaffolds[degreeOfNeighbourhood].begin(); it!= belongingScaffolds[degreeOfNeighbourhood].end(); it++){
//                (*it)->actualize_scaffold(m_orientationRef,indexInBelongingScaffold[degreeOfNeighbourhood].at(*indexIterator));
//                indexIterator++;
//            }
//
//            //Actualize the scaffolds of destinationSite
//            indexIterator=destinationSite.indexInBelongingScaffold[degreeOfNeighbourhood].begin();
//            for (it = destinationSite.belongingScaffolds[degreeOfNeighbourhood].begin(); it!= destinationSite.belongingScaffolds[degreeOfNeighbourhood].end(); it++){
//                (*it)->actualize_scaffold(destinationSite.m_orientationRef,destinationSite.indexInBelongingScaffold[degreeOfNeighbourhood].at(*indexIterator));
//                    indexIterator++;
//            }
//    }
//
    
//    destinationSite.setFullSite(m_orientationRef);
//    this->setFullSite(orientationCopy);
}
     

pair<bool,float> Site::attemptSwapTo(Site &destinationSite, float Temperature, int new_state){
   
    pair<bool, float> successAndDeltaE = {true, 0} ;//the boolean represents the success of swap, and the int is delta_energy
    if (m_orientationRef>0 && destinationSite.isEmpty()){
        
        int old_state = m_orientationRef ;
         
        float old_energy = getSiteEnergy() + destinationSite.getSiteEnergy();
        
//        cout << "Old energy swap: "<<old_energy;
        this->putParticleIn(destinationSite, new_state);
        float new_energy = getSiteEnergy() + destinationSite.getSiteEnergy();
//        cout << "  New energy swap: "<<new_energy<<endl;
        successAndDeltaE.second=new_energy-old_energy;
         if (successAndDeltaE.second>0){
             
             float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

             if (exp (- static_cast <float> (successAndDeltaE.second/Temperature)) < r){
                 destinationSite.putParticleIn(*this, old_state);
                 successAndDeltaE.first=false;
                 successAndDeltaE.second=0;
             }
         }
        
    }
    else{
        cout<< "The initial site was empty or the destination site was not. It can not be swapped" <<endl;
    }
        
//    cout << " success DE "<< sucessAndDeltaE.first << " "<< sucessAndDeltaE.second<< endl;
    return successAndDeltaE;
}

pair <bool, float> Site::attemptFlip(int new_state, float Temperature){

    pair <bool, float> successAndDeltaE {true, 0} ;
    if (!(this->isEmpty())){
        
        float old_orientationRef = m_orientationRef;
        float old_energy = this->getSiteEnergy();
        
        string OldScaffolds = "";
        vector <Scaffold *> :: iterator it;
        for (it = belongingScaffolds[0].begin(); it!= belongingScaffolds[0].end(); it++){
            OldScaffolds += to_string((*it)->getStructure());
        }
        
//        cout << "Old energy flip: "<<old_energy;
        
//        this->fillWithNewRandomValue(numberParticleType, m_orientationRef);
        this->setFullSite(new_state);
        float new_energy = this->getSiteEnergy();
//        cout << " New energy flip: "<<new_energy<<endl;
        successAndDeltaE.second = new_energy - old_energy;
        if (successAndDeltaE.second>0){
            
            float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
            if (exp (- static_cast <float> (successAndDeltaE.second/Temperature)) < r){
                this->setFullSite(old_orientationRef);
                successAndDeltaE.second=0;
                successAndDeltaE.first =false;
            }
        }
 
    }
        else{
        cout<< "The site was empty. It can not be flipped" <<endl;
        }
//    cout << " success DE "<< successAndDeltaE.first << " "<< successAndDeltaE.second<< endl;
    
        return successAndDeltaE;
}


vector <Site*> Site::commonNeigbhours(Site &otherSite, int degreeOfNeighbourhood) const {
    vector <Site *> :: const_iterator it_this ;
    vector <Site *> :: const_iterator it_other ;
    vector <Site *> result;
    for (it_this = m_neighbourList[degreeOfNeighbourhood].begin(); it_this != m_neighbourList[degreeOfNeighbourhood].end() ; ++it_this){
        for (it_other = (otherSite.m_neighbourList[degreeOfNeighbourhood]).begin(); it_other !=(otherSite.m_neighbourList[degreeOfNeighbourhood]).end() ;
            it_other++){
            if ((*it_this)==(*it_other)){
                result.push_back(*it_this);
            }
        }
}
    return result;
}

vector <Scaffold*> Site::commonScaffold(Site &otherSite, int degreeOfNeighbourhood) const{
    vector <Scaffold *> :: const_iterator it_this ;
    vector <Scaffold *> :: const_iterator it_other ;
    vector <Scaffold *> result;
    for (it_this = belongingScaffolds[degreeOfNeighbourhood].begin(); it_this != belongingScaffolds[degreeOfNeighbourhood].end() ; ++it_this){
        for (it_other = (otherSite.belongingScaffolds)[degreeOfNeighbourhood].begin(); it_other !=(otherSite.belongingScaffolds)[degreeOfNeighbourhood].end() ;
            it_other++){
            if ((*it_this)==(*it_other)){
                result.push_back(*it_this);
            }
        }
}
    return result;
}

bool Site::isEmpty() const{
    return m_orientationRef==0;
}

bool Site::equals(Site siteToCompare){
    return (this == &siteToCompare);
}

bool operator==(Site &site1, Site &site2){
    return site1.equals(site2);
}

bool operator!=(Site &site1, Site &site2){
    return !site1.equals(site2);
}
/*
vector <int> Site::getNeighbourOrientationRef () const {
    long N_neighbour = m_neighbourList.size();
    vector <int> NeighbourOrientation (N_neighbour);
    for (int i=0 ; i<N_neighbour ; i++)
    {
        NeighbourOrientation[i]=(*(m_neighbourList[i])).m_orientationRef ;
    }
    return NeighbourOrientation ;
}*/


bool Site::isNeighbourOf(Site& site2, int degreeOfNeighbourhood) const {
    //vector <Site*> &listNeighbour (site2.m_neighbourList);
    bool answ ;
    if (site2.m_neighbourList[degreeOfNeighbourhood].size()==0){
        
        answ = false;
    }
    else{
        answ=(std::find(site2.m_neighbourList[degreeOfNeighbourhood].begin(), site2.m_neighbourList[degreeOfNeighbourhood].end(), this) != site2.m_neighbourList[degreeOfNeighbourhood].end()) ;
    }
    return answ;
}

bool Site::isSecondNeighbour(Site& site2, int degreeOfNeighbourhood) const{
    vector <Site *> :: const_iterator i ;
    bool answ (false); //default is false, whenever find a second neighbour relation, return true.
    for (i = (site2.m_neighbourList)[degreeOfNeighbourhood].begin(); i != (site2.m_neighbourList)[degreeOfNeighbourhood].end() ; ++i){
        if (this->isNeighbourOf((**i),degreeOfNeighbourhood) ){
            answ=true;
            break;
        }
    }
    return answ ;
}


