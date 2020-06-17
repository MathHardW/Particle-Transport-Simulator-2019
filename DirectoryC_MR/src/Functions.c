//===========================IMPORTAÇÕES================================
#include <math.h>
#include <omp.h>
#include <cblas.h>
//======================================================================
#include "FuncoesMatematicas.c"
#include "Alocar.c"
//======================================================================

//===========================STRUCT'S===================================
typedef struct { //DADOS INICIAIS DO PROBLEMA
		int numRegioes;  		    // define o campo número de regiões
		int numMateriais;   		   // define o campo número de zonas materiais

		double *TamanhoRegiao; //define o campo tamanho por região
		int *ParticaoRegiao;    	 // define o campo partição por região
		double *FonteRegiao;    // define o campo fonte por região
		int *ZonaRegiao;       // define o campo zona por região

		double *ST;    // define o campo sigma total por zona
		double *SS;   // define o campo sigma de espalhamento por zona

		double CriterioParada;    // define o campo erro aceitável para o critério de parada
		int OrdemQuadratura;   // define o campo ordem de quadratura

		double ContornoEsquerda;   // define o campo condição de contorno a esquerda
		double ContornoDireita;  // define o campo condição de contorno a esquerda

		int Periodicidade;      // define o campo Periodicidade para exibição dos resultados
} DadosEntrada;// Define o nome do novo tipo criado

typedef struct { // DADOS RESULTADOS DO PROBLEMA
	 int J;

	 AutosVV   *autosVV;
	  double  **MInInversa; // duas dimensões | R e OQ*OQ |
	  double  **psiPart;
      double ***psi;
      double ***MR;
      double *fluxoEscalar;
      
      double ***matrizF;
      double  **psiIn;
} Resultados;

typedef struct { // Cria uma STRUCT para armazenar os dados de entrada do problema
    double **psi;   // define o campo fluxo angular
    double *fi;   // define o campo fluxo escalar
    double *x;   // define o campo posição do domínio para o fluxo escalar
    double *taxaAbsorNodo;// define o campo taxa de absorção por nodo
    double *taxaAbsorRegiao;// define o campo taxa de absorção por região
    double fugaDir;// define o campo fuga pelo contorno direito
    double fugaEsq;// define o campo fuga pelo contorno esquerdo
    double tempoVar; //define o campo tempo de execução das varreduras
    int numInter;// define o campo número de iterações
    int dimFi; // define o campo dimensão do vetor de fluxo escalar
} DiamodDifference;
//======================================================================

//===========================CABEÇALHOS=================================
void ConveterViaSplit(double *ponteiroDouble, char *ponteiroChar);
void *PrepararMatriz(double ***matrizMR, DadosEntrada DadosIniciais);
void matrizMR(Resultados resultados, DadosEntrada DadosIniciais);
void *SolucaoParticular(double **psiPart, DadosEntrada DadosIniciais);

void Varredura(Resultados resultados, DadosEntrada DadosIniciais);
//void _getMult(int n, int j, double *map, int OQ, double ***MR, double **psiPart, double **psi, int dir);
double getMult(int n, int j, int *map, int OQ, Resultados resultados);

DiamodDifference getDD(DadosEntrada d);
void salvarTXTResultado(Resultados r);
//======================================================================

//=============================FUNÇÕES==================================

//FUNÇÃO SPLIT/CONVERSOR================================================
void ConveterViaSplit(double *ponteiroDouble, char *ponteiroChar){
	char temp[100]; int contNum=0, contTemp=0;

	for(int i=0; i<=(sizeof(ponteiroChar)+sizeof(ponteiroChar)); i++){
		if(*(ponteiroChar+i) == ' '){
			*(ponteiroDouble + contNum) = atof(temp);
			//printf("TEMP = %s e Double = %lf\n\n", temp, *(ponteiroDouble + contNum));
			contNum++;
			contTemp=0;
		}else{
			temp[contTemp] = *(ponteiroChar+i);
			contTemp++;
		}
	}
}
//======================================================================

//PREPARAÇÃO DA MATRIZ MR===============================================
void *PrepararMatriz(double ***matriz, DadosEntrada DadosIniciais){
	PLegendre resultado = GaussLegendreAbsPes(DadosIniciais.OrdemQuadratura);
	double *abscissa = resultado.mi;
	double *peso = resultado.w;

	for(int z = 0; z < DadosIniciais.numMateriais; z++){
		for(int j = 0; j < DadosIniciais.OrdemQuadratura; j++){
			for(int k = 0; k < DadosIniciais.OrdemQuadratura; k++){
				if(j==k){
					matriz[z][j][k] = (DadosIniciais.ST[z] - 0.5*DadosIniciais.SS[z]*peso[k])/abscissa[j];
				}else{
					matriz[z][j][k] = (- 0.5*DadosIniciais.SS[z]*peso[k])/abscissa[j];
				}
			}
		}
	}

	free(abscissa);
	free(peso);
}
//======================================================================

//MATRIZ MR ============================================================
void matrizMR(Resultados resultados, DadosEntrada DadosIniciais){
	   int      OQ = DadosIniciais.OrdemQuadratura;
	   int      NR = DadosIniciais.numRegioes;

	  double  *MIn = (double *) malloc(OQ * OQ * sizeof(double));
	  double *MOut = (double *) malloc(OQ * OQ * sizeof(double));

	int cont;
	for (int r = 0; r < NR; r++) {

        double h   = DadosIniciais.TamanhoRegiao[r] / DadosIniciais.ParticaoRegiao[r];
		   int ZMR = DadosIniciais.ZonaRegiao[r];

		cont=0;
		for(int l=0; l<OQ; l++){
			for(int c=0; c<OQ; c++){

				if(l<OQ/2){
					if(resultados.autosVV[ZMR].AutoValor[c]>0){
						MIn[cont] = resultados.autosVV[ZMR].AutoVetor[l][c];
						MOut[cont] = resultados.autosVV[ZMR].AutoVetor[l][c] * exp(-h * resultados.autosVV[ZMR].AutoValor[c]);
					}else{
						MIn[cont] = resultados.autosVV[ZMR].AutoVetor[l][c]  * exp(h * resultados.autosVV[ZMR].AutoValor[c]);
						MOut[cont] = resultados.autosVV[ZMR].AutoVetor[l][c];
					}
				}else{
					if(resultados.autosVV[ZMR].AutoValor[c]>0){
						MIn[cont]  = resultados.autosVV[ZMR].AutoVetor[l][c] * exp(-h * resultados.autosVV[ZMR].AutoValor[c]);
						MOut[cont] = resultados.autosVV[ZMR].AutoVetor[l][c];
					}else{
						MIn[cont]  = resultados.autosVV[ZMR].AutoVetor[l][c];
						MOut[cont] = resultados.autosVV[ZMR].AutoVetor[l][c] * exp(h * resultados.autosVV[ZMR].AutoValor[c]);
					}
				}

				cont++;
			}

		}

		inverse(MIn, OQ);

		double 	*C = (double *) calloc(OQ * OQ , sizeof(double));
		multiplicarMatrizes(OQ, MOut, MIn, C);


		//TRANSFORMANDO O VETOR RESULTADO NA MATRIZ [l][c]
		cont = 0;
		for(int l=0; l<OQ; l++){
			for(int c=0; c<OQ; c++){
				resultados.MR[r][l][c] = C[cont];
				resultados.MInInversa[r][cont] = MIn[cont];
				cont++;
			}
		}

		free(C);
	}

}

//CONSTRUÇÃO DA SOLUÇÃO PARTICULAR======================================
void *SolucaoParticular(double **psiPart, DadosEntrada DadosIniciais){
	int NR = DadosIniciais.numRegioes;
	int ZM = DadosIniciais.numMateriais;
	int OQ = DadosIniciais.OrdemQuadratura;

	PLegendre resultado = GaussLegendreAbsPes(OQ);
	double *matrizLambdaTemp = (double *) malloc(OQ * OQ * sizeof(double));

	double *abscissa = resultado.mi;
	double *peso = resultado.w;

	int OQquad = OQ*OQ , c =0, l=0;;

	for(int r = 0; r < NR; r++){
		int ZMR = DadosIniciais.ZonaRegiao[r];

		//CONSTRUIR LAMBDA [forma eficiente]
		l=0; c=0;
		for(int i = 0; i < OQquad; i++){
			if(c==OQ){c=0; l++;}
			if(l==c){
				matrizLambdaTemp[i] = (DadosIniciais.ST[ZMR] - 0.5 * DadosIniciais.SS[ZMR] * peso[c]);
			}else{
				matrizLambdaTemp[i] = (-0.5*DadosIniciais.SS[ZMR]  *peso[c]);
			}
			c++;
		}

		//INVERSA DA lambda (MATRIZ DE DUAS DIMENSÕES)
		inverse(matrizLambdaTemp, OQ);

		//multiplicar pela fonte
		double tmpSoma=0; l=0; c=0;
		for(int i = 0; i < OQquad; i++){
			if(c==OQ){
				psiPart[r][l] = DadosIniciais.FonteRegiao[r] * tmpSoma;
				tmpSoma = 0;
				l++; c=0;
			}
			tmpSoma += matrizLambdaTemp[i];
			c++;
		}
		psiPart[r][l] = DadosIniciais.FonteRegiao[r] * tmpSoma;

	}

	free(abscissa);
	free(peso);
}
//======================================================================

void Varredura(Resultados resultados, DadosEntrada DadosIniciais){
	int NR = DadosIniciais.numRegioes;
	int ZM = DadosIniciais.numMateriais;
	int OQ = DadosIniciais.OrdemQuadratura;
	
	/* CÁLCULO DA MATRIZ psiIN --------------------------------------*/{
		for(int r=0; r < NR; r++){
			for(int l = 0; l < OQ; l++){
				for(int c = 0; c < OQ/2; c++){
					resultados.psiIn[r][c] = resultados.MR[r][l][c];
					resultados.psiIn[r][c+OQ/2] = resultados.MR[r][l][OQ-c-1];
				}
			}
		}		
	}
	
	/* CÁLCULO DA MATRIZ F ----------------------------------------- */{
		//LAÇO DE REPETIÇÃO PARA O PRIMEIRO INDICE - DE REGIÕES ------------
		for(int r = 0; r<NR; r++){
			//ESTA É A MATRIZ F --------------------------------------------
			double 	*F = (double *) calloc(OQ * OQ , sizeof(double));
			//ESTA É A ALOCAÇÃO DA MATRIZ INTERNA DA R ---------------------
			double 	*matrizInterna = (double *) calloc(OQ * OQ , sizeof(double));
			double 	*psiPart = (double *) calloc(OQ * OQ , sizeof(double));
			//CONTADOR PARA INCREMENTO - MATRIZ DE DUAS DIMENSÕES VIRANDO VETOR
			int cont = 0;
			
			//LAÇO PARA AS LINHAS ------------------------------------------
			printf("\nMATRIZ I - R (%i)\n", r);
			for(int l = 0; l < OQ; l++){
				//LAÇO PARA AS COLUNAS -------------------------------------
				for(int c = 0; c < OQ; c++){
					//CONDIÇÃO PARA DIMINUIR A MATRIZ IDENTIDADE PELA MR ---
					if(l==c){
						matrizInterna[cont] = 1 - resultados.MR[r][l][c];
					}else{
						matrizInterna[cont] = resultados.MR[r][l][c];
					}
					printf("%e   ", matrizInterna[cont]);
					//CONVERTENDO PSIPART DE MATRIZ PARA VETOR ------------- pistPart[r][l]
					psiPart[cont] = resultados.psiPart[l][c];
					cont++;
				}
				printf("\n");
			}
			
			//MULTIPLICANDO A MATRIZ INTERNA DA R PELO PSIPART -------------
			//E GUARGANDO O VALOR NA MATRIZ F ------------------------------
			multiplicarMatrizes(OQ, matrizInterna, psiPart, F);
			
			
			
			
			//aqui embolo ...
			
			double 	*matrizR = (double *) calloc(OQ * OQ , sizeof(double));
			double 	*psiIn = (double *) calloc(OQ * OQ , sizeof(double));
			double 	*RMultIn = (double *) calloc(OQ * OQ , sizeof(double));
			
			//TRANSFORMANDO O F[l][c] no matrizF[r][l][c]
			cont = 0;
			for(int l=0; l<OQ; l++){
				for(int c=0; c<OQ; c++){
					resultados.matrizF[r][l][c] = F[cont];
					
					//convertendo para multiplicação entre R e PsiIN
					matrizR[cont] = resultados.MR[r][l][c];
					psiIn[cont] = resultados.psiIn[r][c];
					
					cont++;
				}
			}
			
			multiplicarMatrizes(OQ, matrizR, psiIn, RMultIn);
			
			cont = 0;
			for(int l=0; l<OQ; l++){
				for(int c=0; c<OQ; c++){
					resultados.psi[r][l][c] = RMultIn[cont] + F[cont];
					cont++;
				}
			}

			free(F); free(matrizInterna);
			free(matrizR); free(psiIn);
			free(RMultIn);
		}
	}
	
}

DiamodDifference getDD(DadosEntrada e){
	DiamodDifference dd;
	/*
	 * Preenche os vetores com os dados do domínio
	 */
	int J = 0;
	for (int i = 0; i < e.numRegioes; i++) J += e.ParticaoRegiao[i];
	int N = e.OrdemQuadratura;
	dd.dimFi = J + 1;
	//NNR = e.NNR;
	//double Stop = e.erro;
	//TipoContorno = e.CCETipo + e.CCDTipo;
	double *St = AlocarDouble1D(J);
	double *Ss = AlocarDouble1D(J);
	double *Q = AlocarDouble1D(J);
	double *h = AlocarDouble1D(J);
	double *S = AlocarDouble1D(J);
	//double *Fi = DoubVet1D(J + 1);
	//double *x = DoubVet1D(J + 1);
	//double *TaxaAbR = DoubVet1D(e.numReg);
	double **psiOld = AlocarDouble2D(e.OrdemQuadratura, J + 1);
	double **psi = AlocarDouble2D(e.OrdemQuadratura, J + 1);
	double tmpErro = 0;
	double erro = 0;
	int k = 0;
	dd.fi = AlocarDouble1D(J + 1);
	dd.x = AlocarDouble1D(J + 1);
	dd.x[0] = 0;
	for (int i = 0; i < e.numRegioes; i++)
	{
		double temp = e.TamanhoRegiao[i] / e.ParticaoRegiao[i];
		for (int j = 0; j < e.ParticaoRegiao[i]; j++)
		{
			h[k] = temp;
			Ss[k] = e.SS[e.ZonaRegiao[i]];
			St[k] = e.ST[e.ZonaRegiao[i]];
			Q[k] = e.FonteRegiao[i];
			dd.x[k + 1] = dd.x[k] + temp;
			k++;
		}
	}
	/*  Estabelecendo a condição de contorno */
	for (int n = 0; n < e.OrdemQuadratura/2; n++)
	{
		if(e.ContornoEsquerda < 0){
			psi[n][0] = 0;
			psiOld[n][0] = psi[n][0];
		}
		else{
			psi[n][0] = e.ContornoEsquerda;
			psiOld[n][0] = psi[n][0];
		}
		if(e.ContornoDireita < 0){
			psi[e.OrdemQuadratura/2 + n][J] = 0;
			psiOld[e.OrdemQuadratura/2 + n][0] = psi[e.OrdemQuadratura/2 + n][0];
		}
		else{
			psi[e.OrdemQuadratura/2 + n][J] = e.ContornoDireita;
			psiOld[e.OrdemQuadratura/2 + n][0] = psi[e.OrdemQuadratura/2 + n][0];
		}
	}
	/* Obtendo as quadraturas de Gauss-Legendre */
	PLegendre pl = GaussLegendreAbsPes(e.OrdemQuadratura);
	double *Mi = pl.mi;
	double *W = pl.w;
	//#############################################
	double start;
	start = omp_get_wtime();
	dd.numInter = 0;
	do{
		dd.numInter++;
		/* Calculando a fonte de espalhamento */
		for (int j = 0; j < J; j++)
		{
			S[j] = 0;
			for (int n = 0; n < N; n++)
				S[j] += 0.25 * Ss[j] * (psi[n][j + 1] + psi[n][j]) * W[n];
		}

		/* Varredura para a direita/esquerda */
		for (int m = 0; m < N / 2; m++)
			for (int j = 0; j < J; j++)
			{
				/* Varredura para a direita */
				psiOld[m][j + 1] = psi[m][j + 1];
				psi[m][j + 1] = ((Mi[m] / h[j] - 0.5 * St[j]) * psi[m][j] + S[j] + Q[j]) / (Mi[m] / h[j] + 0.5 * St[j]);
				/* Varredura para a esquerda */
				psiOld[m + N / 2][J - 1 - j] = psi[m + N / 2][J - 1 - j];
				psi[m + N / 2][J - 1 - j] = ((Mi[m] / h[J - 1 - j] - 0.5 * St[J - 1 - j]) * psi[m + N / 2][J - j] + S[J - 1 - j] + Q[J - 1 - j]) / (Mi[m] / h[J - 1 - j] + 0.5 * St[J - 1 - j]);
			}

		/*  Estabelecendo a condição de contorno */
		for (int n = 0; n < e.OrdemQuadratura/2; n++)
		{
			if(e.ContornoEsquerda < 0){
				psiOld[n][0] = psi[n][0];
				psi[n][0] = psi[e.OrdemQuadratura/2 + n][ 0];
			}
			if(e.ContornoDireita < 0){
				psiOld[e.OrdemQuadratura/2 + n][0] = psi[n][0];
				psi[e.OrdemQuadratura/2 + n][J] = psi[n][J];
			}
		}

		/* Calculando o erro para o critério de parada*/
		tmpErro = 0;
		erro = 0;
		for (int m = 0; m < N; m++)
			for (int j = 0; j < J; j++){
				tmpErro = Abs(psi[m][j] - psiOld[m][j])/psi[m][j];
				if(tmpErro > erro) erro = tmpErro;
			}
	}while(erro > e.CriterioParada);
	dd.tempoVar = omp_get_wtime() - start;

	/* Calcula o fluxo escalar */
	for (int j = 0; j <= J; j++){
		dd.fi[j] = 0;
		for (int m = 0; m < N; m++)
			dd.fi[j] += psi[m][j] * W[m];
	}

	return dd;
}

void salvarTXTResultado(Resultados r){
	FILE *file;
	file = fopen("Result1.txt", "w");

	fprintf(file, "//FLUXO ESCALAR \n");
	for(int i=0; i 	<= r.J; i++){
		fprintf(file, "%lf \n", r.fluxoEscalar[i]);
	}

	fclose(file);
}
