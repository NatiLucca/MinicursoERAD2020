#include <iostream>
#include <omp.h>
#include <vector>
#include <random>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <ostream>
#include"../include/pso.hpp"
#include"../include/benchmarking.hpp"

using namespace std;

extern string nameArq;

PSO::PSO(int population,double boundMin,double boundMax,double C1,double C2,int iterations, double inertia, int dim, string tf){
    this->population = population;
    this->dimension = dim;
    this->rangemin = boundMin;
    this->rangemax = boundMax;
    for(int i=0;i<this->dimension;i++)
      this->global_best_pos.push_back(0);
    this->c1 = C1;
    this->c2 = C2;
    this->steps = iterations;
    this->inertia = inertia;
    this->global_best = 0;
    this->typeFunction = tf;
}
double PSO::getRangeMax(){
  return this->rangemax;
}
double PSO::getRangeMin(){
    return this->rangemin;
}
double PSO::getGlobalBest(){
      return this->global_best;
}
void PSO::initParticles(){
  /* Inicializa todas as particulas em posições aleatorias e velocidade 0 */
  ofstream archiveParticle;
  archiveParticle.open("Inputs/valuesParticleInit" + nameArq + ".txt", ios::ate);
        for (int i =0; i<population; i++){
            Particle P(this->dimension, this->rangemin, this->rangemax, archiveParticle);
            this->particles.push_back(P);
        }
    archiveParticle.close();
}

void PSO::initValuesParticles(){
  ifstream arqParticles;
  arqParticles.open("Inputs/valuesParticleInit" + nameArq + ".txt");
  for (int i =0; i<population; i++){
        Particle P(this->dimension, this->rangemin, this->rangemax, arqParticles);
        this->particles.push_back(P);
  }
    arqParticles.close();
}

/* Caso o custo for menor que o armazenado pelo particula
  atualiza o custo
  atualiza a posição da particula que gerou o custo */
void PSO::melhor_local_global(double cost_Particle, int i, int j){
    if(cost_Particle <= particles[j].best_fitness || i==0){
        particles[j].best_fitness = cost_Particle;
        particles[j].particle_best_position = particles[j].position;
    }
	/*custo global*/
    if(cost_Particle <= this->global_best ){
        this->global_best = cost_Particle;
        this->global_best_pos = particles[j].position;
    }
}

void PSO::calc_velocity_position(double rand1, double rand2, int k, int j){
      double local_diference = (- particles[j].position[k] + particles[j].particle_best_position[k]);
      double global_diference = (- particles[j].position[k] + this->global_best_pos[k]);

      particles[j].velocity[k] = (particles[j].velocity[k]* this->inertia) + rand1 * this->c1 * local_diference + rand2 * this->c2 * global_diference;

      if((particles[j].position[k] + particles[j].velocity[k]) < this->getRangeMax() ){
            particles[j].position[k] += particles[j].velocity[k];
      }else{
            particles[j].position[k] = this->getRangeMax();
      }
}

void PSO::run(double (*costFunction)(vector<double> )){
      ifstream archiveRand;
      string name = "Inputs/valuesPSO" + nameArq + ".txt";
      archiveRand.open(name);
      long double rand1 = 0;
      long double rand2 = 0;
      long double cost_Particle = 0;
      archiveRand >> rand1;
      archiveRand >>rand2;
      archiveRand.close();

     initValuesParticles();

     /* inicializa o melhor global */
     this->global_best = costFunction(this->particles[0].position);

      for (int i=0; i<steps; i++){
             for(int j=0; j<population; j++){
                    cost_Particle = costFunction(particles[j].position);

                    melhor_local_global(cost_Particle, i , j);

                    for(int k = 0; k < dimension ; k++){
                            calc_velocity_position(rand1, rand2, k, j);
                    }
            }/* end for */
      }/* end for steps */
}

void PSO::runParallel(double (*costFunction)(vector<double> )){
      ifstream archiveRand;
      string name = "Inputs/valuesPSO" + nameArq + ".txt";
      archiveRand.open(name);
      long double rand1 = 0;
      long double rand2 = 0;
      long double cost_Particle = 0;
      archiveRand >> rand1;
      archiveRand >>rand2;
      archiveRand.close();

     initValuesParticles();

     /* inicializa o melhor global */
     this->global_best = costFunction(this->particles[0].position);

      for (int i=0; i<steps; i++){
             for(int j=0; j<population; j++){
                    cost_Particle = costFunction(particles[j].position);

                    melhor_local_global(cost_Particle, i , j);

                    for(int k = 0; k < dimension ; k++){
                            calc_velocity_position(rand1, rand2, k, j);
                    }
            }/* end for */
      }/* end for steps */
}

void PSO::runDetails(double (*costFunction)(vector<double> )){

  ofstream archive;
  string name = "Inputs/valuesPSO" + nameArq+ ".txt";
  archive.open(name);

  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<> dis(0.0,1.0);
  long double rand1 = dis(gen);
  long double rand2 = dis(gen);
  long double cost_Particle = 0;
  archive << rand1 << endl << rand2 << endl;
  archive.close();

  initParticles();

}
