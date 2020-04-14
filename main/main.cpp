
#include<iostream>
#include<stdlib.h>
#include <string>
#include <time.h>
#include"../include/particle.hpp"
#include"../include/pso.hpp"
#include"../include/benchmarking.hpp"
#include <sys/time.h>

using namespace std;

string nameArq;

int help(){
  cout << "Error!\n" << "Report\n  population min max c1 c2 steps function inertia dimension [min|max] \n or \n population min max c1 c2 steps function inertia dimension [min|max] nome flag" << endl;
  exit(0);
}
int main(int argc, char** argv){

  int teste =0;
  if(argc == 12){
      nameArq = 1;
      teste =  2;
  }else if(argc != 13 ){
    help();
  }else{
    nameArq =  (argv[11]);
    teste =  atoi(argv[12]);
  }

  int population = atoi(argv[1]);
  double boundMin = atof(argv[2]);
  double boundMax = atof(argv[3]);
  int c1 = atoi(argv[4]);
  int c2 = atoi(argv[5]);
  int steps = atoi(argv[6]);
  string function = (argv[7]);
  double inertia =  atof(argv[8]);
  int dim = atoi(argv[9]);
  string typeFunction = (argv[10]);

  if(typeFunction != "min" && typeFunction != "max"){
    help();
  }

  clock_t start_t, end_t;
  double total_t;
  start_t = clock();

  PSO pso(population, boundMin, boundMax, c1, c2, steps, inertia, dim, typeFunction);

   /* GERA ARQUIVO DE INPUTS*/
  if( teste == 0){
            if (function == "griewank"){
                          pso.runDetails(griewank);
            }else if (function == "alpine"){
                          pso.runDetails(alpine);
            }else if (function == "booth"){
                          pso.runDetails(booth);
             }else if (function == "easom"){
                           pso.runDetails(easom);
            }else if (function == "rosenbrock"){
                          pso.runDetails(rosenbrock);
            }else if (function == "rastringin"){
                          pso.runDetails(rastringin);
            }else if (function == "sphere"){
                          pso.runDetails(sphere);
            }else{
                cout << "\nUnidentified Function" << endl;
          }
  }else if (teste == 1){ /* Versão Sequencial */
          if (function == "griewank"){
                        pso.run(griewank);
          }else if (function == "alpine"){
                        pso.run(alpine);
          }else if (function == "booth"){
                        pso.run(booth);
           }else if (function == "easom"){
                         pso.run(easom);
          }else if (function == "rosenbrock"){
                        pso.run(rosenbrock);
          }else if (function == "rastringin"){
                        pso.run(rastringin);
          }else if (function == "sphere"){
                        pso.run(sphere);
          }else{
              cout << "\nUnidentified Function" << endl;
          }

  }else if (teste == 2){ /* Versão Paralela */
                if (function == "griewank"){
                              pso.runParallel(griewank2);
                }else if (function == "alpine"){
                              pso.runParallel(alpine2);
                }else if (function == "booth"){
                              pso.runParallel(booth2);
                 }else if (function == "easom"){
                               pso.runParallel(easom2);
                }else if (function == "rosenbrock"){
                              pso.runParallel(rosenbrock2);
                }else if (function == "rastringin"){
                              pso.runParallel(rastringin2);
                }else if (function == "sphere"){
                              pso.runParallel(sphere2);
                }else{
                    cout << "\nUnidentified Function" << endl;
                }
  }

  end_t = clock();
  total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;

  /* Saida formato CSV .
  Melhor Solução , tempo de execução */
  cout << pso.getGlobalBest() << " ; " << total_t << endl;

  return 0;
}
