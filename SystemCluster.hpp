//
//  SystemCluster.hpp
//  FrustratedSystem2.0
//
//  Created by lara koehler on 02/12/2020.
//  Copyright © 2020 lara koehler. All rights reserved.
//

#ifndef SystemCluster_hpp
#define SystemCluster_hpp

#include <stdio.h>
#include <vector>
#include <map>
#include <set>
//#endif /* SystemCluster_hpp */

#include "System.hpp"

class SystemCluster {
    
    public :
    SystemCluster(System mySystem);
    
    
    void perform_cluster_measurements(int NClusterMeasurements, int MC_steps_between_measurements, float Temperature_statistics) ; // Does everything in the class
    
    void do_extra_MC_steps(int MC_steps_between_measurements, float Temperature_statistics) ; // Do some more MC steps at the end of the annealing to be sure that the state that is plotted is the one where the energy is minimum
    
    private :
    void reset_cluster_measurements() ; // all the arguments about the cluster measurements are reset to zero

    
    void create_maps(); // Create the two maps that relate each sites to an indicator that will be dealt with in this class
    void create_full_empty_lists() ; // Fills the two toBeTreated lists
    
   
    
    void add_to_current_cluster(int &siteID);
    //    this function will remove the considered site from toBeTreated, add its neighbours in waitingSites, and put it in clusters[clusterInExploration]
    
    void add_to_current_hole(int &siteID);
    //    this function will remove the considered site from toBeTreatedEmpty, add its neighbours in waitingSites, and put it in holes[holeInExploration]
    
    void search_clusters();
    // builds the cluster list from scratch
    
    void save_clusters();
    
    void search_holes();
    // for each cluster of the clusters list, find the emtpy sides that are holes in this cluster

    void save_holes();
    
    void add_a_surface_pair(int &fullSiteID, int &emptySiteID);
    // this pair of full-empty sites that are neighbor is a surface, the lists emptyInSurface, fullInSurface, indicesOfAppartenanceHole and indicesOfAppartenanceCluster will be actualized by adding the information about this new pair (if it was not there already)
    
    void actualize_a_pair(int &emptySiteID);
    // when an empty surface site is found while actualizing the holes references, we update the list of indicesOfAppartenanceHole to add to which hole this empty site belongs, because it was not known when the list was initially filled in the function add_a_surface_pair
    
    void relate_hole_to_cluster();
    
    void summarize_cluster_statistics();
    
    void save_cluster_information() ;
    
    
    
    
    

    
    int cluster_measurement_index ; // the iteration of the measurement of cluster information for a given system, this will be put in the name of the file for cluster information
    
    System mySystem ;
    
    std::map<Site*, int> mapSiteToId ;
    std::map<int, Site*> mapIdToSite ;
    
    std::set<int> toBeTreated ; // the set of indicators of full sites that will be progressively emptied while clusters is filled
    std::set<int> toBeTreatedEmpty ; // the set indicators of empty sites that will be progressively emptied while clusters is filled
    
    
    int clusterInExploration ; // the index of the cluster that is currently being explored
    int holeInExploration ; // the index of the hole that is currently being explored
    
    std::set <int> waitingSites ; // the indicator of sites that will belong to the cluster that is currently being explored
    
    std::vector< std::vector <int>> clusters ; // the final list of clusters
    
    
    std::vector< std::vector <int>> holes ; // list i of this vector is the list of the indicator of the emtpy sites that are holes of cluster i
    int Nparticles;
    
    std::vector <int> emptyInSurface ; //list of indicator of empty sites that are at the surface, a site can appear several times in this list if it has several contacts with an full site
    std::vector <int> fullInSurface ;//list of indicator of full sites that are at the surface, a site can appear several times in this list if it has several contacts with an empty site
    std::vector <int> indicesOfAppartenanceHole ; //indicesOfAppartenanceHole[i] is the index of the hole to which emptyInSurface[i] belongs
    std::vector <int> indicesOfAppartenanceCluster ; //indicesOfAppartenanceCluster[i] is the index of the cluster to which fullInSurface[i] belongs

    std::set <std::pair<int, int>> treatedSurfacePairs ; //set in which we put the indicator of the two sites (full, empty) have already been added to the emptyInSurface and fullInSurface lists

    
    std::vector <int> clusterOfHole ; // clusterOfHole[i] = j means that the hole of index i is within the cluster of index j. If clusterOfHole[i] = -1 it means that hole i is actually the vacuum. The list is initialized with -2 everywhere.
    
    std::vector <int> vacuumID ; //the indices of the holes that are the vacuum
    // they can be unconnected because of the diagonal contacts in some geometries
    
    bool vacuumFound ; //true if the index of the hole that is the solvant (or vacuum) has already been identified 
    
    std::vector <int> Nfull ; // Nfull[k] = nfull means the k-th cluster has nfull particles
    std::vector <int> Nvacancies ; // Nvacancies[k] = nvacancy means the k-th cluster has nvacancy topologically separated holes
    std::vector <int> SizeVacancies ; // SizeVacancies[k] = svacancy means the k-th cluster has svacancy empty sites within the cluster
    std::vector <int> InnerSurface ; // InnerSurface[k] = isurf means the k-th cluster has isurf bonds toward a vacancy
    std::vector <int> OuterSurface ; // OuterSurface[k] = osurf means the k-th cluster has osurf bonds toward a the solvent
    
    
};

#endif /* SystemCluster_hpp */
