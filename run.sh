
#!/bin/bash
# chmod +x run.sh

if [ ! -d "Resultados" ]; then
	mkdir Resultados
else
	rm Resultados/*.csv
fi

if [ ! -d "Inputs" ]; then
	mkdir Inputs
else
	rm Inputs/*.txt
fi
#-fopt-info-vec-optimized
g++ -std=c++11 -o psos main/main.cpp main/benchmarking.cpp main/benchmarkingP.cpp include/benchmarking.hpp src/particle.cpp src/pso.cpp include/particle.hpp include/pso.hpp

g++ -O3 -ftree-vectorize  -mavx -std=c++11 -o pso main/main.cpp main/benchmarking.cpp main/benchmarkingP.cpp include/benchmarking.hpp src/particle.cpp src/pso.cpp include/particle.hpp include/pso.hpp -fopenmp

pop=10
inertia=0.7
typeF="min"
c1=2
c2=2

## BENCHMARK Griewank##
for i in $(seq 1 30);
do
	./pso $pop -600 600 $c1 $c2 5000 griewank $inertia 50 $typeF $i 0 >> temp.csv
done

for i in $(seq 1 30);
do
	./psos $pop -600 600 $c1 $c2 5000 griewank $inertia 50 $typeF $i 1 >> Resultados/GriewankS.csv
done

for i in $(seq 1 30);
do
	./pso $pop -600 600 $c1 $c2 5000 griewank $inertia 50 $typeF $i 2 >> Resultados/GriewankP.csv
done

## BENCHMARK Alpine##
for i in $(seq 1 30);
do
 	./pso $pop 0 10 $c1 $c2 3000 alpine $inertia 10 $typeF $i 0 >> temp.csv
 done

for i in $(seq 1 30);
do
	./psos $pop 0 10 $c1 $c2 3000 alpine $inertia 10 $typeF $i 1 >> Resultados/AlpineS.csv
done

for i in $(seq 1 30);
do
 	./pso $pop 0 10 $c1 $c2 3000 alpine $inertia 10 $typeF $i 1 >> Resultados/AlpineP.csv
done

## BENCHMARK Booth##
for i in $(seq 1 30);
do
 	./pso $pop -10 10 $c1 $c2 1000 booth $inertia 2 $typeF $i 0 >> temp.csv
done

for i in $(seq 1 30);
do
 	./psos $pop -10 10 $c1 $c2 1000 booth $inertia 2 $typeF $i 1 >> Resultados/BoothS.csv
done

for i in $(seq 1 30);
do
 	./pso  $pop -10 10 $c1 $c2 1000 booth $inertia 2 $typeF $i 2 >> Resultados/BoothP.csv
done

## BENCHMARK Easom ##
for i in $(seq 1 30);
do
 	./pso $pop -100 100 $c1 $c2 1000 easom $inertia 2 $typeF $i 0 >> temp.csv
done

for i in $(seq 1 30);
do
 	./psos $pop -100 100 $c1 $c2 1000 easom $inertia 2 $typeF $i 1 >> Resultados/EasomS.csv
done

for i in $(seq 1 30);
do
 	./pso  $pop -100 100 $c1 $c2 1000 easom $inertia 2 $typeF $i 2 >> Resultados/EasomP.csv
done

## BENCHMARK Rosenbrock##
for i in $(seq 1 30);
do
 	./pso $pop -5 10 $c1 $c2 5000 rosenbrock $inertia 50 $typeF $i 0 >> temp.csv
done

for i in $(seq 1 30);
do
 	./psos $pop -5 10 $c1 $c2 5000 rosenbrock $inertia 50 $typeF $i 1 >> Resultados/RosenbrockS.csv
done

for i in $(seq 1 30);
do
 	./pso $pop -5 10 $c1 $c2 5000 rosenbrock $inertia 50 $typeF $i 2 >> Resultados/RosenbrockP.csv
done

## BENCHMARK Rastringin##
for i in $(seq 1 30);
do
  	./pso $pop -5.12 5.12 $c1 $c2 5000 rastringin $inertia 50  $typeF $i 0 >> temp.csv
done

for i in $(seq 1 30);
do
  	./psos $pop -5.12 5.12 $c1 $c2 5000 rastringin $inertia 50  $typeF $i 1 >> Resultados/RastringinS.csv
done

for i in $(seq 1 30);
do
  	./pso $pop -5.12 5.12 $c1 $c2 5000 rastringin $inertia 50  $typeF $i 2 >> Resultados/RastringinP.csv
done

## BENCHMARK Sphere ##
for i in $(seq 1 30);
do
 	./pso $pop -5.12 5.12 $c1 $c2 3000 sphere $inertia 6 $typeF $i 0 >> temp.csv
done

for i in $(seq 1 30);
do
 	./psos $pop -5.12 5.12 $c1 $c2 3000 sphere $inertia 6 $typeF $i 1 >> Resultados/SphereS.csv
done

for i in $(seq 1 30);
do
 	./pso $pop -5.12 5.12 $c1 $c2 3000 sphere $inertia 6 $typeF $i 2 >> Resultados/SphereP.csv
done

rm temp.csv
rm pso
rm psos

exit
