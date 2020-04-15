
#include<iostream>
#include<stdlib.h>
#include <string>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <istream>
#include <ostream>
#include <fstream>

using namespace std;

double media( double s[], int n ){
    double sum = 0.0;
    int i = 0;
    for( i = 0; i < n; i++ )
        sum += s[i];

    return sum / n;
}

double variancia( double s[], int n ){
    double sum = 0.0, dev = 0.0, med = media( s, n );
    int i = 0;
    for( i = 0; i < n; i++ ){
        dev = s[i] - med;
        sum += (dev * dev);
    }
    return sum / n;
}

double desvio_padrao( double s[], int n  ){
    double v = variancia( s, n );
    return sqrt( v );
}


int main(int argc, char** argv){

  string nome = argv[1];
  char x[100] = "";
   double TempS[30],  TempP[30];
   double SolS[30], SolP[30];
   char *p;

   ifstream arq("Resultados/" + nome + "S.csv");
  if(!arq){
    cout << "Não foi possivel abrir o arquivo" << endl;
    return 0;
  }
   int i = 0;
   while(!arq.eof()) {
            arq >> x ; // Solução
            SolS[i] = atof(x);
            arq >> x ; /* Ponto e Virgula ; */
            arq >> x ; // Tempo
            TempS[i] = atof(x);
            i++;
  }
  arq.close();

ifstream arq2("Resultados/" + nome + "P.csv");
 if(!arq2){
   cout << "Não foi possivel abrir o arquivo" << endl;
   return 0;
 }
 i = 0;
while(!arq2.eof()) {
           arq2 >> x ; // Solução
           SolP[i] = atof(x);
           arq2 >> x ; /* Ponto e Virgula ; */
           arq2 >> x ; // Tempo
           TempP[i] = atof(x);
           i++;
 }
 arq2.close();

double  medTempS, medTempP, desvS, desvP, sped;

medTempS = media(TempS, 30); // Media do Tempo sequencial
medTempP = media(TempP, 30); // Media do Tempo paralelo
desvS = desvio_padrao(TempS, 30); // Desvio Padrão do Tempo sequencial
desvP = desvio_padrao(TempP, 30); // Desvio Padrão  do Tempo paralelo
sped = ((medTempS - medTempP) / medTempP ) *100 ;
cout << "Sequencial    -   Paralelo" << endl;
cout << "Média do Tempo de Execução: " << medTempS << "   -   " << medTempP << endl;
cout << "Taxa de Desempenho: " << sped << "%" << endl;

  return 0;
}
