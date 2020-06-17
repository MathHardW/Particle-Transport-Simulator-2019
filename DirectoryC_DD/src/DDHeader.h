#ifndef DDHEADER_H_
#define DDHEADER_H_
//#include <conio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#define MAX_THREADS 12

/*
 * Definindo uma Struct para os dados de entrada
 */
typedef struct // Cria uma STRUCT para armazenar os dados de entrada do problema
{
    int numReg; // define o campo número de regiões
    int numZon; // define o campo número de zonas materiais

    double *tamReg; // define o campo tamanho por região
    int *partReg; // define o campo partição por região
    double *fontReg; // define o campo fonte por região
    int *zonReg; // define o campo zona por região

    double *stZon; // define o campo sigma total por zona
    double *ssZon; // define o campo sigma de espalhamento por zona

    double erro; // define o campo erro aceitável para o critério de parada
    int ordQuad; // define o campo ordem de quadratura
    double cce; // define o campo condição de contorno a esquerda
    double ccd; // define o campo condição de contorno a esquerda
    int p; // define o campo Periodicidade para exibição dos resultados

    int J;
} DadosEntrada; // Define o nome do novo tipo criado

/*
 * Definindo uma Struct e uma função para as quadraturas de Gauss-Legandre
 */
typedef struct // Cria uma STRUCT para armazenar os dados de entrada do problema
{
    double *mi; // define o campo direção mi
    double *w; // define o peso w
} PLegendre; // Define o nome do novo tipo criado
PLegendre GaussLegendreAbsPes(int n);
double Abs(double x);

/*
 * Manipulação de arquivos
 */
typedef struct // Cria uma STRUCT para armazenar os dados de entrada do problema
{
    char **linha; // define o campo condição de contorno a esquerda
    int dim; // define o campo condição de contorno a esquerda
} DadosLinha; // Define o nome do novo tipo criado
DadosEntrada Ler(char nomeArquivo[]);
DadosLinha LinhaSplit(char str[], char split[]);

/*
 * Manipulação de vetores e matrizes de números inteiros
 */
int *IntVet1D(int dim1);
int **IntVet2D(int dim1, int dim2);
void LiberarInt1D(int *vet1D);
void LiberarInt2D(int **vet2D, int dim1);

/*
 * Manipulação de vetores e matrizes de números reais
 */
double *DoubVet1D(int dim1);
double **DoubVet2D(int dim1, int dim2);
double ***DoubVet3D(int dim1, int dim2, int dim3);
void LiberarDoub1D(double *vet1D);
void LiberarDoub2D(double **vet2D, int dim1);
void LiberarDoub3D(double ***vet2D, int dim1, int dim2);

/*
 * Definindo uma Struct e as funções para a execução do Método Diamod Difference
 */
typedef struct // Cria uma STRUCT para armazenar os dados de entrada do problema
{
    double **psi; // define o campo fluxo angular
    double *fi; // define o campo fluxo escalar
    double *x; // define o campo posição do domínio para o fluxo escalar
    double *taxaAbsorNodo; // define o campo taxa de absorção por nodo
    double *taxaAbsorRegiao; // define o campo taxa de absorção por região
    double fugaDir; // define o campo fuga pelo contorno direito
    double fugaEsq; // define o campo fuga pelo contorno esquerdo
    double tempoVar; //define o campo tempo de execução das varreduras
    int numInter; // define o campo número de iterações
    int dimFi; // define o campo dimensão do vetor de fluxo escalar
} DiamodDifference;
DiamodDifference getDD(DadosEntrada d);
DiamodDifference getDDParallel(DadosEntrada d);

void Escrever(DadosEntrada dados, DiamodDifference DD);

#endif /* DDHEADER_H_ */

/*
printf("%d\n",dados.numReg);
printf("%d\n",dados.numZon);
for(int i = 0; i < dados.numReg; i++){printf("%.2f\t",dados.tamReg[i]);}printf("\n");
for(int i = 0; i < dados.numReg; i++){printf("%d\t",dados.partReg[i]);}printf("\n");
for(int i = 0; i < dados.numReg; i++){printf("%.2f\t",dados.fontReg[i]);}printf("\n");
for(int i = 0; i < dados.numReg; i++){printf("%d\t",dados.zonReg[i]);}printf("\n");
for(int i = 0; i < dados.numZon; i++){printf("%.2f\t",dados.stZon[i]);}printf("\n");
for(int i = 0; i < dados.numZon; i++){printf("%.2f\t",dados.ssZon[i]);}printf("\n");
printf("%.2f\n",dados.erro);
printf("%d\n",dados.ordQuad);
printf("%.2f\n",dados.cce);
printf("%.2f\n",dados.ccd);
 */

//printf("Diamond Diferene Serial\nForam realizada %d iterações, com tempo de varreduras igual a %.5f segundos\n\n", DD.numInter, DD.tempoVar);
//printf("Diamond Diferene Paralelo\nForam realizada %d iterações, com tempo de varreduras igual a %.5f segundos\n\n\n", DDParallel.numInter, DDParallel.tempoVar);

//printf("Eficiencia de %.3f por cento",100*(DD.tempoVar - DDParallel.tempoVar)/DD.tempoVar);

//int n = 8;
//PLegendre p = GaussLegendreAbsPes(n);
//for(int i = 0; i < n; i++){printf("Mi[%d] = %.10f\tW[%d] = %.10f\n",i+1,p.mi[i],i+1,p.w[i]);}printf("\n");
