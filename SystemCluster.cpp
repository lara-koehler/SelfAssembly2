//
//  SystemCluster.cpp
//  FrustratedSystem2.0
//
//  Created by lara koehler on 02/12/2020.
//  Copyright © 2020 lara koehler. All rights reserved.
//

#include "SystemCluster.hpp"
#include "LogDuration.hpp"


#include "System.hpp"


#include <time.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <set>

using namespace std;
    

SystemCluster::SystemCluster(System arg_mySystem) :
cluster_measurement_index(0),
mySystem(arg_mySystem),
mapSiteToId(map<Site*, int>()),
mapIdToSite(map<int, Site*>()),
toBeTreated(set<int>()),
toBeTreatedEmpty(set<int>()),
clusterInExploration(0),
holeInExploration(0),
waitingSites(set <int>()),
clusters(vector<vector<int>> (0)),
holes(vector<vector<int>>(0)),
Nparticles(0),
emptyInSurface(vector<int>()),
fullInSurface(vector<int>()),
indicesOfAppartenanceHole(vector<int>()),
indicesOfAppartenanceCluster(vector<int>()),
treatedSurfacePairs(set<pair<int, int>>()),
clusterOfHole(vector<int>()),
vacuumID(vector<int>()),
vacuumFound(false),
Nfull(vector<int>()),
Nvacancies(vector<int>()),
SizeVacancies(vector<int>()),
InnerSurface(vector<int>()),
OuterSurface(vector<int>())
{

    
//    cout << "Initialisation of SystemCluster"<<endl;
    
    
    create_maps();
    
    Nparticles = (int)toBeTreated.size();
    
//    cout << "end of initialisation"<<endl;
    
}

void SystemCluster::reset_cluster_measurements(){
    toBeTreated = set<int>();
    clusterInExploration = 0 ;
    toBeTreatedEmpty = set<int>() ;
    holeInExploration = 0;
    holes = vector<vector<int>>(0) ;
    waitingSites = set <int>() ;
    clusters = vector<vector<int>> (0);
    emptyInSurface = vector<int>();
    fullInSurface= vector<int>() ;
    indicesOfAppartenanceHole = vector<int>();
    indicesOfAppartenanceCluster = vector<int>();
    treatedSurfacePairs = set<pair<int, int>>() ;
    clusterOfHole = vector<int>() ;
    Nfull = vector<int>();
    Nvacancies = vector<int>() ;
    SizeVacancies = vector<int>() ;
    InnerSurface = vector<int>() ;
    OuterSurface = vector<int>() ;
    vacuumID = vector<int>() ;
    vacuumFound = false ;
}


void SystemCluster::perform_cluster_measurements(int NClusterMeasurements, int MC_steps_between_measurements, float Temperature_statistics){
    LOG_DURATION("Cluster measurements ("+to_string(NClusterMeasurements)+" measures)");

    
    for (int t=0 ; t<NClusterMeasurements ; t++){
        
        create_full_empty_lists();
//        cout<<"Lists created, ";
        search_clusters();
//        cout << "Clusters measured, ";
        search_holes();
//        cout << "Holes measured, ";
        
        // The two following functions will provide informations about the hole and the clusters adress but they are only useful for debug
//        save_clusters();
//        save_holes();
        
        
        relate_hole_to_cluster();
//        cout<<"Holes related to clusters, ";
        summarize_cluster_statistics();
//        cout<< "Statistics summarized, ";
        save_cluster_information();
//        cout<<"Information saved, "<<endl;
        
        reset_cluster_measurements();
//        cout<<"Reset done"<<endl;
        cluster_measurement_index ++ ;
        if (t!=NClusterMeasurements-1){
            do_extra_MC_steps(MC_steps_between_measurements, Temperature_statistics);
        }
    }
    
}

void SystemCluster::create_maps(){
// fills mapSiteToId and mapIdToSite that relate each site to its indicator
    for (int i=0; i<mySystem.get_Nsites() ; i++){
        Site* aSite = mySystem.get_site(i);
        mapSiteToId[aSite] = i;
        mapIdToSite[i] = aSite ;
        
        }
}

void SystemCluster::create_full_empty_lists(){
    for (int i=0; i<mySystem.get_Nsites() ; i++){
    Site* aSite = mySystem.get_site(i);
    if (aSite->isEmpty()) {
        // indicators of empty sites are added to the set toBeTreatedEmpty
        int emptySiteID = i ;
        toBeTreatedEmpty.insert(emptySiteID) ;
    }
    else {
        // indicators of full sites are added to the set toBeTreated
        int fullSiteID = i ;
        toBeTreated.insert(fullSiteID) ;
    }
    }
}

    

void SystemCluster::do_extra_MC_steps(int MC_steps_between_measurements, float Temperature_statistics){
    // run a few monte-carlo steps to reach ground state (state with the minimal energy measured)
    
    
    double Ecurrent = mySystem.get_Etot() ;
    double Emin = mySystem.get_Emin() ;
    int nMCstep = 0 ;
    int Nmes = MC_steps_between_measurements;
    
    // do a few steps to change a bit the configuration
    for (int i=0 ; i<Nmes; i++){
        mySystem.MC_step(Temperature_statistics, false);
    }
    
    // do enough steps to recover a state of minimal energy
    while (nMCstep<Nmes and Ecurrent>Emin){
        mySystem.MC_step(Temperature_statistics, false);
        nMCstep++ ;
    }
    
    // print wether the minimum of energy was reached
//    if (nMCstep == Nmes){
//        cout << "Emin "<< Emin <<" was not reached, "<< mySystem.get_Etot() <<" instead "<<endl;
//    }
//    else{
//    cout << "Emin "<<Emin<<" was reached in "<<nMCstep<<" out of "<<    Nmes<< endl ;
//    }
    
}


void SystemCluster::add_to_current_cluster(int &siteID){
    //cout <<"the site we consider "<<endl;
    //aSite.printAttributes(3);
    
    if (not mapIdToSite[siteID]->isEmpty()){
        // the site is full
        for (auto& neighbour : mapIdToSite[siteID]->getNeighbours(0)){
            
            if (not neighbour->isEmpty()){
                // the neighbour is full
                int neigbourID = mapSiteToId[neighbour] ;
                if (toBeTreated.find(neigbourID)!=toBeTreated.end()){
                    //check if the neighbour has not be treated yet, else add it to the same cluster later
                    waitingSites.insert(neigbourID);
                }
            }
            else{
                // the neighbour is empty
                int neigbourID = mapSiteToId[neighbour] ;
                // we have found a surface pair, we will deal with it in add_a_surface_pair
                add_a_surface_pair(siteID, neigbourID);
                }
        }
        int newSiteID = siteID ;
        toBeTreated.erase(newSiteID);
        clusters[clusterInExploration].push_back(newSiteID);
    }
   
}

void SystemCluster::add_a_surface_pair(int &fullSiteID, int &emptySiteID){
//    A surface pair has been found, if it was not in the list already, the reference of the two sites and the index of the cluster to which the full site belongs are pushed back in the corresponding lists
    pair<int, int > pairUnderStudy;
    pairUnderStudy.first = fullSiteID ;
    pairUnderStudy.second = emptySiteID ;
    
    if (treatedSurfacePairs.find(pairUnderStudy)==treatedSurfacePairs.end()) {
        // if it has not been treated already
        treatedSurfacePairs.insert(pairUnderStudy);
        emptyInSurface.push_back(emptySiteID);
        fullInSurface.push_back(fullSiteID);
        indicesOfAppartenanceCluster.push_back(clusterInExploration);
        indicesOfAppartenanceHole.push_back(holeInExploration);
        
    }
}

void SystemCluster::actualize_a_pair(int &surfaceEmptySiteID){
    for (int index=0; index<emptyInSurface.size(); index++){
//        every time the emptysite considered appears in the list of sites that are in a surface, we actualize the value of the hole to which it belongs
        if (surfaceEmptySiteID==emptyInSurface[index]){
//            aSurfaceEmptySite is in the pair of index i, and now we can actualize to which hole it belongs to
            indicesOfAppartenanceHole[index] = holeInExploration;
        }
    }
}



void SystemCluster::add_to_current_hole(int &siteID){
    if (mapIdToSite[siteID]->isEmpty()){
        bool isOnSurface = false ;
        for (auto& neighbour : mapIdToSite[siteID]->getNeighbours(0)){
            if ( neighbour->isEmpty()){
                int neigbourID = mapSiteToId[neighbour] ;
                if (toBeTreatedEmpty.find(neigbourID)!=toBeTreatedEmpty.end()){
                //check if the neighbour has not be treated yet, else :
                    waitingSites.insert(neigbourID);
                }
            }
            else{
//                if we find at least one non-empty neighbor, it means we are dealing with a surface site, and we need to actualize the value of the hole to which it belongs in indicesOfAppartenanceHole. This will be done with the function actualize_pair
                isOnSurface = true ;
            }
        }
        if (isOnSurface){
            actualize_a_pair(siteID);
        }
        
//        int newSiteID = siteID ;
        toBeTreatedEmpty.erase(siteID);
        holes[holeInExploration].push_back(siteID);
    }
}



void SystemCluster::search_clusters(){
    int count = 0;
    while (not toBeTreated.empty()){
        count+=1;
        clusters.push_back(vector <int> (0));
        // we add a new list to clusters, which correspond to the current cluster
        
        int seedSiteID = *(toBeTreated.begin()) ;
        //seedSite will be the seed of the new cluster
        add_to_current_cluster(seedSiteID);
        //seedSite is removed from toBeTreated, its neighbours are added to waitingSites,
        // and it is added to clusters
        int count2=0;
        while (not waitingSites.empty()){
            count2++;
            int aSiteID = *(waitingSites.begin());
            add_to_current_cluster(aSiteID);
            int aSiteIDtoErase = aSiteID;
            waitingSites.erase(aSiteIDtoErase);
        }
        clusterInExploration+=1;
        //cout << "size of toBeTreated in search_cluster "<<toBeTreated.size()<<endl;
        //cout << "clusters size in search_cluster "<<clusters.size()<<endl;
    }
        
    }

void SystemCluster::search_holes(){
    int count = 0;
    while (not toBeTreatedEmpty.empty()){
        count+=1;
        holes.push_back(vector <int> (0));
        // we add a new list to clusters, which correspond to the current cluster
        int seedSiteID = *(toBeTreatedEmpty.begin());
        //seedSite will be the seed of the new cluster
        add_to_current_hole(seedSiteID);
        //seedSite is removed from toBeTreated, its neighbours are added to waitingSites,
        // and it is added to clusters
        int count2=0;
        while (not waitingSites.empty()){
            count2++;
            int aSiteID = *(waitingSites.begin());
            add_to_current_hole(aSiteID);
            waitingSites.erase(aSiteID);
        }
        holeInExploration+=1;
        //cout << "size of toBeTreated in search_cluster "<<toBeTreated.size()<<endl;
        //cout << "clusters size in search_cluster "<<clusters.size()<<endl;
    }
}


void SystemCluster::relate_hole_to_cluster(){
    
    // Initialize the clusterOfHole vector
    for (int h=0 ; h<holeInExploration; h++){
        int holeID_temp (-2);
        clusterOfHole.push_back(holeID_temp);
    }
    
    
    if (holeInExploration==1){
        //There is only one hole in the system, it is necessarily the vacuum
        clusterOfHole[0] = -1;
        vacuumID.push_back(0) ;
        vacuumFound = true ;
    }
    
    else{
    
    for (int clusterID=0 ; clusterID<indicesOfAppartenanceHole.size() ; clusterID++){
        int holeID = indicesOfAppartenanceHole[clusterID];

        if (clusterOfHole[holeID] == -2){
            // nothing has been done with this hole yet
            clusterOfHole[holeID] = indicesOfAppartenanceCluster[clusterID]; ;
            // we relate the hole to its corresponding cluster
        }
        else {
            //if (clusterOfHole[holeID]!= -1 and clusterOfHole[holeID]!=indicesOfAppartenanceCluster[clusterID]){
              // this assumes that there is only one vacuum. However with square geometry, vacuum can be isolated by corner contacts, but this is still the vacuum and not a vacancy
            
                if (clusterOfHole[holeID]!=indicesOfAppartenanceCluster[clusterID]){
                //we found a vacuum, it's a hole connected to two different clusters
                clusterOfHole[holeID] = -1;
                vacuumID.push_back(holeID) ;
                vacuumFound = true ;
            }
        }
    }
        if (not vacuumFound){
            // we are in the situation of one single cluster that has holes
            int maxSize = 0 ;
            for (int holeID=0 ; holeID<holes.size() ; holeID++){
                // we go through all the holes to find the bigger one, which is the vacuum
                if (holes[holeID].size() > maxSize ){
                    maxSize = (int) holes[holeID].size() ;
                    vacuumID.push_back(holeID) ;
                    vacuumFound = true ;
                }
            }
            if (vacuumFound){
                clusterOfHole[vacuumID[0]] = -1;
            }
            else{
                cout<< " PROBLEM, vacuum was not found"<<endl;
            }
        }
    }
//    for (int k=0 ; k<clusterOfHole.size(); k++){
//        cout << "hole "<<k<<" in cluster "<<clusterOfHole[k]<<endl;
//    }
}



void SystemCluster::summarize_cluster_statistics(){
    
//    if (clusters.size() == 0){
//        cout << endl<<"THERE IS A PROBLEM IN CLUSTER COMPUTATION, NUMBER OF CLUSTERS IS ZERO"<<endl;
//    }
//    if (holes.size() == 0){
//        cout << endl<<"THERE IS A PROBLEM IN CLUSTER COMPUTATION, NUMBER OF HOLES IS ZERO"<<endl;
//    }
//
    //Initializing the lists
    for (int nc=0 ;  nc<clusters.size() ; nc++){
        int nf_temp = 0 , nv_temp = 0, sv_temp = 0, is_temp = 0, os_temp = 0 ; // initialization of the values in the list (to avoid "double free or corruption (out)" error
        Nfull.push_back(nf_temp) ;
        Nvacancies.push_back(nv_temp) ;
        SizeVacancies.push_back(sv_temp) ;
        InnerSurface.push_back(is_temp) ;
        OuterSurface.push_back(os_temp) ;
    }
    
    // Filling Nfull
    for (int nc=0 ; nc<clusters.size() ; nc++){
        Nfull[nc] = (int)clusters[nc].size();
    }
    
    //Filling Nvacancies
    for (int nh=0 ; nh<holes.size() ; nh++){
        int clusterItBelongsTo = clusterOfHole[nh] ;
        if ( clusterItBelongsTo >= 0){
            // the vacuum (that belongs to cluster -1) is not taken into account
            Nvacancies[clusterItBelongsTo] +=1 ;
            SizeVacancies[clusterItBelongsTo] += holes[nh].size();
        }
//        else {
//            cout << "Hole number "<<nh<< " is assigned to cluster "<<clusterItBelongsTo<<endl;
//        }
        
    }
    
    //Filling InnerSurface and OuterSurface
    for (int link=0 ; link<indicesOfAppartenanceCluster.size(); link++){
        int clusterOfTheLink = indicesOfAppartenanceCluster[link];
        int holeOfTheLink = indicesOfAppartenanceHole[link];
        
        
        //if (holeOfTheLink==vacuumID){
        if (std::find(vacuumID.begin(), vacuumID.end(), holeOfTheLink)!=vacuumID.end()){
            OuterSurface[clusterOfTheLink]+=1 ;
        }
        
        else if (clusterOfHole[holeOfTheLink]==clusterOfTheLink){
            InnerSurface[clusterOfTheLink]+=1 ;
        }
        else{
            cout<<"Problem in summarize_cluster_statistices, hole "<<holeOfTheLink<<" should be in cluster "<<clusterOfHole[holeOfTheLink]<<" and not in "<<clusterOfTheLink <<endl ;
        }
        
        
    }
    
//    for (int nc=0 ;  nc<clusters.size() ; nc++){
//        cout <<"cluster "<<nc<<" has "<<Nfull[nc]<<" full sites, "<<Nvacancies[nc]<<" vacancies, "<<InnerSurface[nc]<< " inner surface, and "<< OuterSurface[nc]<<" outer surface"<<endl;
//    }
   
}



void SystemCluster::save_clusters(){
    string  directory1;
    directory1 = mySystem.get_directoryResults()+"_clusters.txt";
    ofstream myStream1 (directory1.c_str());
    
    vector <vector<int>> :: iterator aCluster;
    vector <int> :: iterator sitesID;
    int countClusters = 0;
    int countSize = 0 ;
    
//    cout <<"Sizes of clusters : ";
    for (aCluster = clusters.begin(); aCluster!=clusters.end(); aCluster++){

        //cout << "Cluster "<<countClusters+1<<" is being explored, ";
        for (sitesID = clusters[countClusters].begin(); sitesID!=clusters[countClusters].end(); sitesID++){
            myStream1<<mapIdToSite[(*sitesID)]<<" ";
            countSize++;
        }
//        cout <<countSize<<" " ;
        countClusters+=1;
        countSize=0;
        myStream1<<endl;
    }
    
    
    cout<< countClusters <<" clusters in this system"<<endl;
    myStream1<<"#END";
}


void SystemCluster::save_holes(){
    string  directory1;
    directory1 = mySystem.get_directoryResults()+"_holes.txt";
    ofstream myStream1 (directory1.c_str());

    vector <vector<int>> :: iterator aHole;
    vector <int> :: iterator sitesID;
    int countHoles = 0;
    int countSize = 0 ;
    cout <<"Sizes of holes : ";
    for (aHole = holes.begin(); aHole!=holes.end(); aHole++){

        //cout << "Cluster "<<countClusters+1<<" is being explored, ";
        for (sitesID = holes[countHoles].begin(); sitesID!=holes[countHoles].end(); sitesID++){
            myStream1<<mapIdToSite[(*sitesID)]<<" ";
            countSize++;
        }
        cout <<countSize<<" " ;
        countHoles+=1;
        countSize=0;
        myStream1<<endl;
    }
    cout <<endl<< countHoles <<" holes in this system"<<endl;
    myStream1<<"#END";
}

void SystemCluster::save_cluster_information(){
    string  directory1;
    directory1 = mySystem.get_directoryResults()+"_Nfull_Nvac_SizeVac_InSurf_OutSurf_"+to_string(cluster_measurement_index)+".txt";
    ofstream myStream1 (directory1.c_str());
    
    int numberOfClusters = (int)clusters.size();
//    cout <<"Sizes of clusters : ";
    
    for (int nc=0 ;  nc<numberOfClusters ; nc++){
//        cout<<Nfull[nc]<<" ";
        myStream1<< to_string(Nfull[nc])<<" "<<to_string(Nvacancies[nc])<<" "<<to_string(SizeVacancies[nc])<<" "<<to_string(InnerSurface[nc])<< " "<< to_string(OuterSurface[nc])<<" "<<mapIdToSite[clusters[nc][0]];
        myStream1<<endl;
    }
    if (cluster_measurement_index==0){
        cout << numberOfClusters <<" clusters in this system"<<endl;
    }
    myStream1<<"#END";
    
}




    

