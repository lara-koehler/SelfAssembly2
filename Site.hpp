//
//  Site.hpp
//  FrustratedSystems
//
//  Created by lara koehler on 22/09/2020.
//  Copyright © 2020 lara koehler. All rights reserved.
//

#ifndef Site_hpp
#define Site_hpp

#include <stdio.h>
#include <string>
#include <vector>
#include "Scaffold.hpp"

class Site {
    public :
    
    Site(); //constructeur
    /*Site(int orientationRef, std::vector<std::vector <Site*>> neighbourList);*/
    
    //Site* getAdress() const;
    
    void setScaffoldAppartenance(Scaffold* s_adress, int indexInS, int degreeOfNeighbourhood);

    
    std::vector <Scaffold*> getBelongingScaffolds(int degreeOfNeighbourhood){
        return belongingScaffolds[degreeOfNeighbourhood];
    };
    
    
    bool setNeighbourRelation(Site& neighbourSite,int degreeOfNeighbourhood); //add the adress of neighbourSite to the attribute list of the current site, and vice versa. Returns true if no neighbour relation existed before
    
    void  printAttributes(int const Nhash) const ; //outputs the adress of a site, its orientation, and the adresses of its neighbour
    
    /*static void printVectorOfSite(std::vector <Site *> listToPrint);*/ //MADE A TEMPLATE METHOD INSTEAD
    
//    void fillWithRandomValue(int numberParticleType);
//    void fillWithNewRandomValue (int numberParticleType, int currentOrientation);
    
    void initializeEmptySite();
    
    void initializeFullRandom ();
    
    void setFullSite (int orientationRef);
    void changeScaffoldOrientation();
    void changeScaffoldStructure();
    
    float getSiteEnergy();
    
    void setOrientation(int orientationValue); // just changes the value of m_orientationRef without changing anything on the scaffolds etc
    void putParticleIn (Site &destinationSite, int new_state); //if site2 is an empty site, put particle of site in destinationSite
        // potentially also changes the state of the particle
    
    std::pair<bool,float> attemptSwapTo(Site &destinationSite, float Temperature, int new_state);
    std::pair<bool,float>  attemptFlip(int numberParticleType, float Temperature);
    
    int getOrientationRef() const;
    int getSiteStateTemporary();
    
    std::vector <int> getNeighbourOrientationRef (int degreeOfNeighbourhood) const;
    
    bool isEmpty() const ;
    
    bool equals(Site siteToCompare) ;
    
    bool isNeighbourOf(Site& site2, int degreeOfNeighbourhood) const ; //check if two sites are neighbours
    /* this is a non symetric relation, testing if this is referenced in the neighbour of site 2.*/
    
    bool isSecondNeighbour(Site& site2, int degreeOfNeighbourhood) const ;
    
    void set_position(int const i, int const  j,int  const  k, int const  N);
    
    int get_hashedPosition(){
        return hashedPosition;
    }
    std::vector<int> get_real_space_position(int const N) const;
    
    
    std::vector <Site*> commonNeigbhours(Site &otherSite,int degreeOfNeighbourhood) const;
    std::vector <Scaffold*> commonScaffold(Site &otherSite, int degreeOfNeighbourhood) const;
    
    std::vector <Site*> getNeighbours (int degree);
    
    
    
    private :
    
    int m_orientationRef;
    std::vector<std::vector <Site*> > m_neighbourList ;// the element at position i on this vector is the list of neighbours of the (i+1)th degree of the site (it can only have one element if we only consider nearest neighbours, two if we consider nearest neighbours and second nearest neighbours, etc..
    int hashedPosition;
    std::vector<std::vector <Scaffold*> > belongingScaffolds ;
    std::vector<std::vector <int> > indexInBelongingScaffold ;
    const int maxDegreeOfNeighbourhood = 2;
    
    
};

#endif /* Site_hpp */
