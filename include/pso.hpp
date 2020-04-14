#ifndef PSO_H
#define PSO_H

#include<math.h>
#include<stdlib.h>
#include<vector>
#include <string>
#include"../include/particle.hpp"

using namespace std;

class PSO{

  private :
    vector<Particle> particles;
    int population;
    int dimension;
    int steps;
    double global_best;
    vector<double> global_best_pos;
    double c1,c2;
    double rangemax;
    double rangemin;
    double inertia;
    string typeFunction;

    void initParticles();
    void initValuesParticles();
    void melhor_local_global(double cost_Particle, int i,int j);
    void calc_velocity_position(double rand1, double rand2, int k, int j);

  public :
    PSO(int population,double boundMin,double boundMax,double C1,double C2,int iterations, double inertia, int dim, string tf);
    void run(double (*costFunction)(vector<double> ));
    void runDetails(double (*costFunction)(vector<double> ));
    double getGlobalBest();
    void runParallel(double (*costFunction)(vector<double> ));
    double getRangeMax();
    double getRangeMin();
};

#endif
