//
//  main.cpp
//  FrustratedSystems
//
//  Created by lara koehler on 22/09/2020.
//  Copyright © 2020 lara koehler. All rights reserved.
//

#include <iostream>
#include "Site.hpp"
#include "Scaffold.hpp"
#include <vector>
#include "System.hpp"

#include "SystemType.hpp"
#include "head.hpp"
#include "LogDuration.hpp"

#include <ctime> // For random
#include <cstdlib> // For random

#include <stdio.h>

#include "json.hpp"
using json = nlohmann::json;

#include "SystemCluster.hpp"
using namespace std ;



int main(int argc, const char * argv[]) {
    
    
    string directoryForResults =  argv[1];
    std::ifstream f(directoryForResults+"input_parameters.json");
    json parameters = json::parse(f);
    //cout <<parameters;
    
    LOG_DURATION("Simulation duration");

    //Create the system
    System sys(parameters) ;
    int firstNeighborDistance = parameters["firstNeighborDistance"];
    int secondNeighborDistance = parameters["secondNeighborDistance"];
    sys.create_system(firstNeighborDistance,secondNeighborDistance);
    
    // Get temperature ramp of the annealing
    int NTemperatures = parameters["Ntemperatures"];
    float Temperature_start = parameters["Temperature_start"];
    float Temperature_end = parameters["Temperature_end"];
    std::string Type_temperature_ramp = parameters["Type_temperature_ramp"];
    std::vector <float> Temperatures;
    if (!Type_temperature_ramp.compare("linear")){
        float a = (Temperature_start-Temperature_end) / (NTemperatures-1) ;
        float b = Temperature_start ;
        for (int i =0 ; i<NTemperatures; i++){
            Temperatures.push_back(-a * i + b);
        }
    }
    else if (!Type_temperature_ramp.compare("logarithmic")){
        float a = - (Temperature_start-Temperature_end) / log(NTemperatures) ;
        float b = Temperature_start ;

        for (int i =0 ; i<NTemperatures; i++){
            Temperatures.push_back((a*log(1+i) +b));
        }
    }
    else{
        cout <<"Type_temperature_ramp should be 'linear' or 'logarithmic', not "<<Type_temperature_ramp<<endl;
    }
    
    // Do the annealing
    sys.anneal_system(Temperatures);
    
    // Do the quenching
    float Temperature_quenching = parameters["Temperature_quenching"];
    int Nsteps_quenching = parameters["Nsteps_quenching"];
    sys.quench_system(Temperature_quenching, Nsteps_quenching);
    
    // Save the equilibration information
    sys.save_equilibration();
    
    // Do the statistics on system composition
    int Naverage = parameters["Naverage"];
    bool saveCovariance = parameters["saveCovariance"] ;
    bool saveExtendedScaffold = parameters["saveExtendedScaffold"] ;
    float Temperature_statistics = parameters["Temperature_statistics"];
    if (parameters["doStatistics"]){
        sys.do_statistics(Temperature_statistics, Naverage, saveCovariance, saveExtendedScaffold);
    }
    
    // Do the cluster measurement
    int NClusterMeasurement = parameters["Ncluster_measurements"];
    int MC_steps_between_measurements = parameters["MC_steps_between_measurements"];
    SystemCluster clusters(sys);
    clusters.perform_cluster_measurements(NClusterMeasurement, MC_steps_between_measurements, Temperature_statistics);
    

    return 0;
}
