#include <iostream>
#include <algorithm>
#include <vector>

using namespace std;

bool compare(pair<float,float> a, pair<float,float> b) {
    
    if (a.first<b.first) return true; else return false;
    
};

float get_integral(vector<pair<float,float>> v_masse){
    float integral=0;
    for (vector<pair<float,float>>::const_iterator i = v_masse.begin(); i != v_masse.end(); ++i)
    {
        pair<float,float> temp = *i;
        integral+=temp.first*temp.second;
    }
    return integral;
};


bool is_ok_window(vector<pair<float,float>> v_masse, int indexdown, float sigmaeff, float soglia){
    int indexup=indexdown;
    
        if (v_masse[indexup].first-v_masse[indexdown]<sigmaeff) return true; else return false;
}

bool sliding_match(std::vector<float> v_masse, float sigmaeff){

    
    int soglia = get_integral(v_masse)*0.683;
    
    
    int indexdown=0;
    int indexup=0;
    
    float intfinestra=0;
    while (intfinestra<soglia && indexup<v_masse.size()){
        indexup++;
        intfinestra+=v_masse[indexup].first*v_masse[indexup].second;
    }
    if (v_masse[indexup].first-v_masse[indexdown]<sigmaeff) return true;
    
    while (indexup<v_masse.size()){
        intfinestra-=v_masse[indexdown].first*v_masse[indexdown].second;
        intfinestra+=v_masse[indexup].first*v_masse[indexup].second;
        indexdown++;
        indexup++;
        
        
    }
    
    
    
    
    
    
    for (indexup = soglia; indexup<nentries; indexup++){
        if (v_masse[indexup]-v_masse[indexdown]<0) {std::cout << "sorting error" << std::endl; return false;}
        if (v_masse[indexup]-v_masse[indexdown]<2*sigmaeff) return true;
        indexdown++;
    }
    
    return false;
    
};

double scanvalues_sigmaeff(const char* filename, int max=30000){
    
    
    std::vector<std::pair<float,float>> v_masse;
    
    ifstream in;
    in.open(filename);
    int n=0;
    while(n<max){
	    float x;
	    float w;
	    in >> x;
	    in >> w;
	    if (!in.good()) break;
	    if (70<x && x<110){
            pair<float,float> temp(x,w);
            v_masse.at(i)=temp;
	    }
	    n++;
    }
    in.close();

  
    
    
    std::sort(v_masse.begin(),v_masse.end(),compare);
    
    
    float sigmaeff = 0;
    while (!sliding_match(v_masse,sigmaeff) && sigmaeff<1000) sigmaeff+=0.01;
    
    return sigmaeff;
    
};

