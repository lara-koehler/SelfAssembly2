//
//  System.hpp
//  FrustratedSystems
//
//  Created by lara koehler on 24/09/2020.
//  Copyright © 2020 lara koehler. All rights reserved.
//

#ifndef System_hpp
#define System_hpp

#include "Scaffold.hpp"
#include "Site.hpp"
#include "SystemType.hpp"
#include "head.hpp"

#include "json.hpp"
using json = nlohmann::json;

#include <vector>
#include <map>
#include <stdio.h>
#include <set>



class System {
    
    public:
    //initialisation methods
    System(json parameters);
    
    void build_lattice(int const Lx, int const Ly, int const Lz);
    
    void build_neighbour_relations(int degreeOfNeighbourhood, int distanceOfRelation, std::vector <Scaffold *> &scaffoldListToFill);
    
    void fill_the_system();
    
    void save_parameters();
    
    std::pair<bool,int>  MC_step(float Temperature, bool saveStat);
    std::pair<bool,int>  MC_step_swap(float Temperature, bool saveStat,
                                      bool isFlipSwap);
    std::pair<bool,int>  MC_step_flip(float Temperature, bool saveStat, bool isMutation);
    
    void create_system(int firstNeighborDistance, int secondNeighborDistance);
    void anneal_system(std::vector<float> Temperatures_ramp);
    void quench_system(float Temperature_quenching, int Nsteps_quenching);
    double get_Etot();
    double get_Emin();
    int get_MC_steps_per_Temperature();
//    float get_final_temperature();
    bool get_saveAnnealing();
    int get_Nsites();
    
    double compute_Etot();
    
    void try_save_ground_state();
    // save the system states if it's the lowest measured energy
    // this function also updates the value of Emin
    
    void save_system_state(int indexOfT, int current_step);
    
    void save_equilibration();
    
    std::vector <int> get_scaffold_list_2sites(Site* s1, Site* s2, int degreeOfNeighborhood);
    std::vector <int> get_scaffold_list_1site(Site* s1, int degreeOfNeighborhood);
    
    std::vector <int>  scaffold_stats(std::vector <Scaffold*> aScaffoldList);
    std::vector <int>  scaffold_stats2(std::vector <Scaffold*> aScaffoldList);
    
    void do_statistics(float Temperature_statistics, int Naverage, bool saveCovariance, bool saveExtendedScaffold);
    
    void test_System_initialization();
    void test_site_swap();
    void test_site_flip();
    void test_whatever();
    
    std::vector<Site*> get_full_sites() const;
    std::vector<Site*> get_empty_sites() const;
    
    Site* get_site(int index) const ;
    
    
    std::string get_directoryResults();
    
    void print_energy_per_particle();
    
    
    
    
    private :
    std::vector <Scaffold*> scaffoldList; //should be inherited from a TopologyInitializer class
    std::vector <Scaffold*> extendedScaffoldList;
    std::vector <Site*> siteList; //should be inherited from a TopologyInitializer class
    std::map <int, int> RealPositionMappedToSiteList; //TO REMOVE, should only belong to the topologyInitializer class
    
    
    std::vector <Site*> fullSites;
    std::vector <Site*> emptySites;
    double Etot ;
    double Emin ; // the minimal energy the system was able to reach during annealing
    
    int Lx; 
    int Ly;
    int Lz;
    int NhashPosition; //TO REMOVE, should only belong to the topologyInitializer class
    int Nsites;
    
    int NparticleTypes ;
    std::vector <int>  Nparticles ;
    int NparticlesTot ;
    int  N_Orientation;
    int  latticeFilename ;
    std::string latticeTypeName;
    int  N_Structures ;
    std::vector <float> LEL ;
    
    float FlipRate ;
    float SwapRate ;
    float FlipSwapRate ;
    float MutationRate ;
    

    
    int Nsteps_per_T;
    int Nmes_per_T;
    int Ntemperatures ;
    
    json parameters ;
    std::string directoryForResults;
    
    bool saveSnapshot;
    bool saveEnergies ;
    bool saveTemperatures ;
    bool saveMCsuccess ;
    
    
    
    bool oneGroundStateSaved ;
    
    std::vector<double> AllEnergies ;
    std::vector<bool> AllMCsuccess ;

    
    std::vector <int> cVector ;
    std::vector <int> dVector ;
    
    
//    std::vector <std::vector<int>> scaffoldsStats ;
//    std::vector <std::vector<int>> extendedScaffoldsStats ;
    
    
    SystemType typeOfThisSystem ;
    
    
    
    
    
};
#endif /* System_hpp */
