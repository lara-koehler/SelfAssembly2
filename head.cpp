//
//  head.cpp
//  FrustratedSystems
//
//  Created by lara koehler on 29/09/2020.
//  Copyright © 2020 lara koehler. All rights reserved.
//

#include "head.hpp"
#include <vector>
#include <iostream>

using namespace std ;

int hash_ints (int const a, int const b, int const c, int const N){
    //cout<<"N"<<N;
    return(a+b*N+c*N*N);
}

int uModulus(int const a, int const N){
    if (N==0){
        return (0);
    }
    if (a>0){
        return(a%N);
    }
    else {
        return((a+N)%N);
    }
}

vector <string> parse(string s, string delimiter){
    size_t pos = 0;
    string token;
    vector <string> result ;
    while ((pos = s.find(delimiter)) != string::npos) {
        token = s.substr(0, pos);
        result.push_back(token);
        //std::cout << token << std::endl;
        s.erase(0, pos + delimiter.length());
    }
    result.push_back(s);
    //std::cout << s << std::endl;
    return result ;
}

int get_int_parameter(string s, string expected_parameter){
    vector <string> keyAndValue = parse(s, " ");
    int res (0);
    if (keyAndValue[0]==expected_parameter){
        res = stoi(keyAndValue[1]);
    }
    else{
        cout<< "Parameter "<<expected_parameter<<" not found in parameter file"<<endl;
    }
    return(res);
}

string get_string_parameter(string s, string expected_parameter){
    vector <string> keyAndValue = parse(s, " ");
    string res;
    if (keyAndValue[0]==expected_parameter){
        res = keyAndValue[1];
    }
    return(res);
}

bool get_bool_parameter(string s, string expected_parameter){
    vector <string> keyAndValue = parse(s, " ");
//    cout <<"get_bool_parameters "<<keyAndValue[0] <<" shows '"<<keyAndValue[1]<<"'"<<endl ;
    bool res=false;
    if (keyAndValue[0]==expected_parameter){
        res = (keyAndValue[1]=="True") ; 
    }
   
    return(res);
}





vector<int> unhash_ints(const int hashedValue, const int key){
    
    vector <int> result ;
    result.push_back(hashedValue % key);
    result.push_back(hashedValue % (key*key) /key);
    result.push_back(hashedValue / (key*key));
    return result;
}
