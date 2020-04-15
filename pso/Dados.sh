#!/bin/bash
# chmod +x dados.sh
if [ ! -d "Performance" ]; then
	mkdir Performance
else
	rm Performance/*.txt
fi

g++ -std=c++11 -o dados speedup.cpp

## BENCHMARK Alpine##
./dados Alpine > Performance/Alpine.txt

## BENCHMARK Booth##
./dados Booth > Performance/Booth.txt

## BENCHMARK Easom ##
./dados Easom > Performance/Easom.txt

## BENCHMARK Griewank##
	./dados Griewank > Performance/Griewank.txt

## BENCHMARK Rosenbrock##
./dados Rosenbrock > Performance/Rosenbrock.txt

## BENCHMARK Rastringin##
./dados Rastringin > Performance/Rastringin.txt

## BENCHMARK Sphere ##
./dados Sphere > Performance/Sphere.txt

rm dados

exit
