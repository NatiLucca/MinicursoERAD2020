
## Compilar

g++ -std=c++11 -o pso main/main.cpp main/benchmarking.cpp main/benchmarkingP.cpp include/benchmarking.hpp src/particlecpp src/pso.cpp include/particle.hpp include/pso.hpp -fopenmp

## Observação
Para a execução paralela é obrigatório o uso de  "-fopenmp", para as demais execuções é opcional.

## Opções de Funções de Teste
<ul>
  <li>alpine</li>
  <li>booth</li>
  <li>easom</li>
  <li>griewank</li>
  <li>rastringin</li>
  <li>rosenbrock</li>
  <li>sphere</li>
</ul>

## Variáveis

population: (int) valor inteiro que corresponde ao número de agentes da população.
boundMin: (double) valor decimal que corresponde ao limite inferior do espaço de busca.
boundMax: (double) valor decimal que corresponde ao limite superior do espaço de busca.
c1: (double) valor decimal que corresponde ao coeficiente social do bando de pássaros.
c2: (double) valor decimal que corresponde ao coeficiente cognitivo do bando de pássaros.
steps: (int) valor inteiro que corresponde ao npumero máximo de iterações do algoritmo.
inertia: (double) valor decimal que corresponde a inércia.
dim: (int) valor inteiro que corresponde a dimensão do enxame.
typeF: (string) campo de texto [min|max] que corresponde a função de minimização
(min) ou maximização (max)

## Objeto PSO
PSO pso(population, boundMin, boundMax, c1, c2, steps, inertia, dim, typeF);

Exemplo:
PSO pso(20, -100, 100, 2, 2, 2000, 0.7, 50, min);

#### Execução Detalhada:
pso.runDetails(function);

Exemplo:
pso.runDetails(griewank);

#### Execução Sequencial:
pso.run(function);

Exemplo:
pso.run(sphere);

#### Execução Paralela
pso.runParallel(function);

Exemplo:
pso.runParallel(rastringin);
