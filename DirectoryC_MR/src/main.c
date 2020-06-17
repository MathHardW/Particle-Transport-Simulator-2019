//===========================IMPORTAÇÕES================================
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//----------------------------------------------------------------------
#include "Functions.c"
//======================================================================

//========================CÓDIGO PRINCIPAL==============================
int main(int argc, char *argv [ ]) {

DadosEntrada DadosIniciais;
Resultados  resultados;
//======================================================================

/*=======================LEITURA DO ARQUIVO==========================*/{
	FILE *file; file = fopen("/home/matheus/Documentos/Projeto_de_Pesquisa/1-Projeto_Oficial/projetoIntegrado/DataFiles/Oficial1", "r"); char frase[100];

	//VERIFICANDO O ARQUIVO DE TEXTO====================================
	if(file == NULL){
		printf("Não foi possível abrir o arquivo!\n");
		getchar();
		exit(0);
	}

	//CAPTURANDO E CONVERTENDO OS DADOS DO ARQUIVO DE TEXTO=============
	int cont=0;
	while(fgets(frase, 1000, file) != NULL){
		if(cont % 2 != 0){
			if(cont == 1){DadosIniciais.numRegioes   = atoi(frase);};
			if(cont == 3){DadosIniciais.numMateriais = atoi(frase);};
			//==========================================================
			if(cont == 5){
				char *ponteiroChar;
				ponteiroChar = (char *) calloc(1, sizeof(frase));
				ponteiroChar = frase;

				double *ponteiroDouble;
				ponteiroDouble = (double *) calloc(DadosIniciais.numRegioes, sizeof(double));

				ConveterViaSplit(ponteiroDouble, ponteiroChar);

				DadosIniciais.TamanhoRegiao = AlocarDouble1D(DadosIniciais.numRegioes);
				for(int i=0; i<DadosIniciais.numRegioes; i++){
					DadosIniciais.TamanhoRegiao[i] = *(ponteiroDouble+i);
				}
			};
			//==========================================================
			if(cont == 7){
				char *ponteiroChar;
				ponteiroChar = (char *) calloc(1, sizeof(frase));
				ponteiroChar = frase;

				double *ponteiroDouble;
				ponteiroDouble = (double *) calloc(DadosIniciais.numRegioes, sizeof(int));

				ConveterViaSplit(ponteiroDouble, ponteiroChar);

				DadosIniciais.ParticaoRegiao = AlocarInteger1D(DadosIniciais.numRegioes);
				for(int i=0; i<DadosIniciais.numRegioes; i++){
					DadosIniciais.ParticaoRegiao[i] = *(ponteiroDouble+i);
				}
			};
			//==========================================================
			if(cont == 9){
				char *ponteiroChar;
				ponteiroChar = (char *) calloc(1, sizeof(frase));
				ponteiroChar = frase;

				double *ponteiroDouble;
				ponteiroDouble = (double *) calloc(DadosIniciais.numRegioes, sizeof(double));

				ConveterViaSplit(ponteiroDouble, ponteiroChar);

				DadosIniciais.FonteRegiao = AlocarDouble1D(DadosIniciais.numRegioes);
				for(int i=0; i<DadosIniciais.numRegioes; i++){
					DadosIniciais.FonteRegiao[i] = *(ponteiroDouble+i);
				}
			};
			//==========================================================
			if(cont == 11){
				char *ponteiroChar;
				ponteiroChar = (char *) calloc(1, sizeof(frase));
				ponteiroChar = frase;

				double *ponteiroDouble;
				ponteiroDouble = (double *) calloc(DadosIniciais.numRegioes, sizeof(double));

				ConveterViaSplit(ponteiroDouble, ponteiroChar);

				DadosIniciais.ZonaRegiao = AlocarInteger1D(DadosIniciais.numRegioes);
				for(int i=0; i<DadosIniciais.numRegioes; i++){
					DadosIniciais.ZonaRegiao[i] = *(ponteiroDouble+i);
				}
			};
			//==========================================================
			if(cont == 13){
				char *ponteiroChar;
				ponteiroChar = (char *) calloc(1, sizeof(frase));
				ponteiroChar = frase;

				double *ponteiroDouble;
				ponteiroDouble = (double *) calloc(DadosIniciais.numMateriais, sizeof(double));

				ConveterViaSplit(ponteiroDouble, ponteiroChar);

				DadosIniciais.ST = AlocarDouble1D(DadosIniciais.numRegioes);
				for(int i=0; i<DadosIniciais.numMateriais; i++){
					DadosIniciais.ST[i] = *(ponteiroDouble+i);
				}
			};
			//==========================================================
			if(cont == 15){
				char *ponteiroChar;
				ponteiroChar = (char *) calloc(1, sizeof(frase));
				ponteiroChar = frase;

				double *ponteiroDouble;
				ponteiroDouble = (double *) calloc(DadosIniciais.numMateriais, sizeof(double));

				ConveterViaSplit(ponteiroDouble, ponteiroChar);

				DadosIniciais.SS = AlocarDouble1D(DadosIniciais.numMateriais);
				for(int i=0; i<DadosIniciais.numMateriais; i++){
					DadosIniciais.SS[i] = *(ponteiroDouble+i);
				}
			};
			//==========================================================
			if(cont == 17){DadosIniciais.CriterioParada   = atoi(frase);};
			if(cont == 19){DadosIniciais.OrdemQuadratura   = atoi(frase);};
			if(cont == 21){DadosIniciais.ContornoEsquerda   = atoi(frase);};
			if(cont == 23){DadosIniciais.ContornoDireita   = atoi(frase);};
		}
		cont++;
	}

	/*
	//IMPRESSÕES DOS RESULTADOS=========================================
	printf("\nNúmero de Regiões    .:  |   %i   | ", DadosIniciais.numRegioes);
	printf("\nNúmero de Materiais  .:  |   %i   | ", DadosIniciais.numMateriais);
	printf("\nOrdem de Quadratura  .:  |   %i   | ", DadosIniciais.OrdemQuadratura);
	printf("\n\nCritério de Parada   .:  | %.3lf | ",DadosIniciais.CriterioParada);
	printf("\nContorno da Esquerda .:  | %.3lf | ", DadosIniciais.ContornoEsquerda);
	printf("\nContorno da Direita  .:  | %.3lf | ", DadosIniciais.ContornoDireita);
	printf("\n\nTamanho por região   .: ");
									for(int i = 0; i < DadosIniciais.numRegioes; i++){
										printf(" | %.3lf | ", DadosIniciais.TamanhoRegiao[i]);
									};
	printf("\nPartição por região  .: ");
									for(int i = 0; i < DadosIniciais.numRegioes; i++){
										printf(" |   %i   | ", DadosIniciais.ParticaoRegiao[i]);
									};
	printf("\nFonte por região     .: ");
									for(int i = 0; i < DadosIniciais.numRegioes; i++){
										printf(" | %.3lf | ", DadosIniciais.FonteRegiao[i]);
									};
	printf("\nMaterial por região  .: ");
									for(int i = 0; i < DadosIniciais.numRegioes; i++){
										printf(" |   %i   | ", DadosIniciais.ZonaRegiao[i]);
									};
	printf("\n\nST por Material      .: ");
									for(int i = 0; i < DadosIniciais.numMateriais; i++){
										printf(" | %.3lf | ", DadosIniciais.ST[i]);
									};
	printf("\nSS por Material      .: ");
									for(int i = 0; i < DadosIniciais.numMateriais; i++){
										printf(" | %.3lf | ", DadosIniciais.SS[i]);
									};
	*/
	fclose(file);

	printf("\n");
}
//======================================================================

int NR = DadosIniciais.numRegioes;
int ZM = DadosIniciais.numMateriais;
int OQ = DadosIniciais.OrdemQuadratura;

resultados.J = 0;
for (int r = 0; r < NR; r++) {
    resultados.J += DadosIniciais.ParticaoRegiao[r];
}

//===========================DECLARAÇÕES================================
double ***matriz = AlocarDouble3D(ZM, OQ, OQ);
PrepararMatriz(matriz, DadosIniciais);
//----------------------------------------------------------------------
resultados.autosVV = CalcAutosVV(matriz,OQ,ZM);
//----------------------------------------------------------------------
resultados.MR = AlocarDouble3D(NR, OQ, OQ);
resultados.MInInversa = AlocarDouble2D(NR, OQ*OQ);
//----------------------------------------------------------------------
matrizMR(resultados, DadosIniciais);
//----------------------------------------------------------------------
resultados.psiPart = AlocarDouble2D(NR, OQ);
//----------------------------------------------------------------------
SolucaoParticular(resultados.psiPart, DadosIniciais);
//----------------------------------------------------------------------
resultados.psi = AlocarDouble3D(NR, OQ, OQ);
resultados.fluxoEscalar = AlocarDouble1D(resultados.J);
//----------------------------------------------------------------------

/*========================IMPRIMIR MATRIZ R============================*/{
printf("\nMATRIZ R [TRÊS DIMENSÕES]===============================\n\n");
	for(int r=0; r<NR; r++){
		for(int l=0; l<OQ; l++){
			for(int c=0; c<OQ; c++){

				if(resultados.MR[r][l][c] < 0){
					printf("| %e | ", resultados.MR[r][l][c]);
				}else{
					printf("|  %e  | ", resultados.MR[r][l][c]);
				}

			}
			printf("\n");
		}
		printf("\n");
	}
}
//======================================================================

resultados.matrizF = AlocarDouble3D(NR, OQ, OQ);
resultados.psiPart = AlocarDouble2D(NR, OQ);
resultados.psiIn = AlocarDouble2D(NR, OQ);
Varredura(resultados, DadosIniciais);
//----------------------------------------------------------------------
DiamodDifference dd = getDD(DadosIniciais);
//----------------------------------------------------------------------
salvarTXTResultado(resultados);
//======================================================================

/*======================IMPRIMIR PSI PARTICULAR======================*/{

printf("\nPSI PARTICULAR=========================================\n\n");
for(int r=0; r<NR; r++){
	for(int l=0; l<OQ; l++){
		printf("|\t %e \t|", resultados.psiPart[r][l]);
	}
	printf("\n\n");
}
printf("\n=======================================================\n\n");
}
//======================================================================

/*========================IMPRIMIR MATRIZ F============================*/{
printf("\nMATRIZ F [TRÊS DIMENSÕES]=============================\n\n");

for(int r=0; r<NR; r++){
	for(int l=0; l<OQ; l++){
		for(int c=0; c<OQ; c++){

			if(resultados.matrizF[r][l][c] < 0){
				printf("| %e | ", resultados.matrizF[r][l][c]);
			}else{
				printf("|  %e  | ", resultados.matrizF[r][l][c]);
			}

		}
		printf("\n");
	}
	printf("\n");
}

printf("\n=======================================================\n\n");
}
//======================================================================

/*======================IMPRIMIR PSI PARTICULAR======================*/{

printf("\nPSI IN=========================================\n\n");
for(int r=0; r<NR; r++){
	for(int l=0; l<OQ; l++){
		printf("|\t %e \t|", resultados.psiIn[r][l]);
	}
	printf("\n\n");
}
printf("\n=======================================================\n\n");
}
//======================================================================


/*========================IMPRIMIR PSI============================*/{
printf("\n PSI - TRÊS DIMENSÕES]===============================\n\n");
	for(int r=0; r<NR; r++){
		for(int l=0; l<OQ; l++){
			for(int c=0; c<OQ; c++){
				if(resultados.psi[r][l][c] < 0){
					printf("| %e | ", resultados.psi[r][l][c]);
				}else{
					printf("|  %e  | ", resultados.psi[r][l][c]);
				}
			}
			printf("\n");
		}
		printf("\n");
	}
}
//======================================================================

/*============================FLUXO ESCALAR=========================={

printf("\nFLUXO ESCALAR MR==========================================\n\n");
for(int j=0; j<=resultados.J; j++){

	if(resultados.fluxoEscalar[j] < 0){
		printf("| %e |", resultados.fluxoEscalar[j]);
	}else{
		printf("| %e  |", resultados.fluxoEscalar[j]);
	}

	printf("\n");
}
printf("\n=======================================================\n\n");

printf("\nFLUXO ESCALAR DD=======================================\n\n");
for(int j=0; j<=resultados.J; j++){

	if(dd.fi[j] < 0){
		printf("| %e |", dd.fi[j]);
	}else{
		printf("| %e  |", dd.fi[j]);
	}

	printf("\n");
}
printf("\n=======================================================\n\n");

}
//======================================================================
*/

return 0;
}

//PARA COMPILAR .: cd "%d"; gcc -o "%e" "%f" -I lib -lopenblas -fopenmp -lm;
