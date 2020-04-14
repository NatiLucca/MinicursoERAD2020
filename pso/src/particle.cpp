#include <iostream>
#include <omp.h>
#include <vector>
#include <random>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <ostream>
#include <cstdlib>
#include"../include/particle.hpp"

using namespace std;


Particle::Particle(){
}

/* Inicializa para cada partícula um número aleatorio e velocidade zero */
Particle::Particle(int dim, double rangemin, double rangemax){
  long double x=0;
      for(int i=0; i<dim; i++){
	/* Produz valores aleatórios de ponto flutuante, distribuídos no intervalo [min, max) */
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> dis(rangemin,rangemax);
        x = dis(gen);
        cout << x << endl;
        this->position.push_back(x);
        this->particle_best_position.push_back(2);
        this->velocity.push_back(0);
      }

}

Particle::Particle(int dim, double rangemin, double rangemax, ostream& arq){
  long double x=0;
      for(int i=0; i<dim; i++){
	/* Produz valores aleatórios de ponto flutuante, distribuídos no intervalo [min, max) */
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> dis(rangemin,rangemax);
        x = dis(gen);
        arq << x << endl;
        this->position.push_back(x);
        this->particle_best_position.push_back(2);
        this->velocity.push_back(0);
      }
}

Particle::Particle(int dim, double rangemin, double rangemax, istream& arq){
  long double x=0;
      for(int i=0; i<dim; i++){
        arq >> x;
        this->position.push_back(x);
        this->particle_best_position.push_back(2);
        this->velocity.push_back(0);
      }
}
