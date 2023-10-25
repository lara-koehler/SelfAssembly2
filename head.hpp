//
//  head.hpp
//  FrustratedSystems
//
//  Created by lara koehler on 29/09/2020.
//  Copyright © 2020 lara koehler. All rights reserved.
//

#ifndef head_hpp
#define head_hpp

#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

int hash_ints (int const a, int const b, int const c, int const N);

int uModulus(int const a, int const N);

std::vector <int> unhash_ints(const int hashedValue, const int key);

std::vector <std::string> parse(std::string s, std::string delimiter) ;

int get_int_parameter(std::string s, std::string expected_parameter) ;
std::string get_string_parameter(std::string s, std::string expected_parameter) ;
bool get_bool_parameter(std::string s, std::string expected_parameter);

template <typename T>
T* remove_at(std::vector<T *>& v, typename std::vector<T>::size_type n)
// Constant-time element eraser for vectors of pointers : first
// switch with last element, then pop_back().
{
    T* to_erase =v.at(n);
    v.at(n)=v.back();
    v.pop_back();
    return to_erase;
}

template <typename T>
void printVectorOfAdress(std::vector <T *> listToPrint){
    //if listToPrint.site()!=0){
    typename std::vector <T *> :: const_iterator i ;
    for (i = listToPrint.begin(); i != listToPrint.end() ; ++i){
        std::cout << *i <<" ";}
    std::cout <<std::endl;
}


template <typename T>
void printVectorOfVariable(std::vector <T> listToPrint){
    //if listToPrint.site()!=0){
    typename std::vector <T> :: const_iterator i ;
    for (i = listToPrint.begin(); i != listToPrint.end() ; ++i){
        std::cout << *i <<" ";}
    std::cout <<std::endl;
}



template <typename T>
void vectorToFile( std::vector<T> vectorToSave, std::string filename){
    typename std::vector <T> :: const_iterator i ;
    std::ofstream myStream (filename.c_str());

        for (i = vectorToSave.begin(); i != vectorToSave.end() ; ++i){
            myStream << (*i) <<" ";
            
    }
    myStream <<"#END";
    myStream.flush();
    myStream.close();
    myStream.clear();
}


template <typename T>
void matrixToFile( std::vector<std::vector <T> > matrixToSave, std::string filename){
    int N1 = (int) matrixToSave.size();
    int N2=0;
    std::ofstream myStream (filename.c_str());
    if (N1!=0) {N2 = (int) matrixToSave[0].size();}
    for (int i=0 ; i<N1 ; i++){
    
        for (int j=0 ; j<N2; j++){
            myStream << matrixToSave[i][j] <<" ";
        }
        myStream<<std::endl;
        
    }
    myStream <<"#END";
    myStream.flush();
    myStream.close();
    myStream.clear();
}




template <typename T>
std::vector <T> VectorFromFile(std::string filename){
    std::vector <T> result ;
    std::fstream myStream(filename.c_str());
    if (myStream) {
    std::string token;
    
    getline(myStream,token);
    std::vector <std::string> string_list (parse(token, " "));
    
  
    int Ns = (int) string_list.size();
    for (int i=0; i<Ns-1;i++){
        T value;
        std::istringstream ss(string_list[i]);
        ss >> value;
        result.push_back(value);
        }
    }
    else {
        std::cout <<"File not found. Filename : " <<filename <<std::endl;
    }
    
    return result ;
}
    


#endif /* head_hpp */
