#ifndef PARTICLEPSO_H
#define PARTICLEPSO_H

#include<math.h>
#include<stdlib.h>
#include <ostream>
#include<vector>

using namespace std;

class Particle{
  private :
    vector<double> velocity; /* velocidade da particula*/
    vector<double>  particle_best_position; /* melhor posição da particula */
    vector<double> position; /* posição x, y da particula */
     double best_fitness;

  public :
    friend class PSO;
    Particle();
    Particle(int dim, double rangemin, double rangemax);
    Particle(int dim, double rangemin, double rangemax, ostream& arq);
    Particle(int dim, double rangemin, double rangemax, istream& arq);

};

#endif
