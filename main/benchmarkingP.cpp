#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <omp.h>
#include "../include/benchmarking.hpp"

using namespace std;

double easom2(vector<double> temp){
    return (-cos(temp[0])) * cos(temp[1]) * exp(-pow(temp[0] - M_PI, 2) - pow(temp[1] - M_PI, 2));
}

double booth2(vector<double> temp){
   return (pow((temp[0]+2*temp[1]-7), 2) + pow((2*temp[0]+temp[1]-5), 2));
}

double sphere2(vector<double> temp){
   double sum=0;
   for(int i=0; i< temp.size(); i++){
	    sum += pow(temp[i], 2);
   }
   return sum;
}

double griewank2(vector<double> temp) {
	double sum=0, prod = 1;
	for (int i = 0; i < temp.size (); i++) {
		sum += pow (temp[i], 2.);
		prod*= cos(temp[i] / sqrt (i+1));
	}
	return (1.0 + (sum / 4000.) - prod);
}

double alpine2(vector<double> temp){
  double sum=0;
  for(int i=0; i< temp.size(); i++){
 	    sum += fabs((temp[i] * sin(temp[i])) + 0.1 * temp[i]);
  }
   return sum;
}

double rastringin2(vector<double> temp){
   double sum=0;
   for(int i=0; i< temp.size(); i++){
 	    sum += pow(temp[i], 2) - 10 * cos(2 * M_PI * temp[i]) + 10;
   }
   return sum;
}

double rosenbrock2(vector<double> temp){
   double sum=0, aux = 0;
   for(int i=0; i< (temp.size()-1); i++){
        aux = (temp[i+1], 2) - pow(temp[i], 2);
 	      sum +=  (100 * pow( aux, 2) + pow((temp[i] - 1), 2));
   }
   return sum;
}
