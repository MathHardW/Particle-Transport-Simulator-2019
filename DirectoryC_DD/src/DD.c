#include "DDHeader.h"

/*
 * Manipulação de arquivos
 */
DadosEntrada Ler(char nomeArquivo[]) {
	FILE *arq;
	char Linha[100];
	char **tmpDados = 0;
	char split[] = " ";
	char *result;
	int l = 0;
	DadosEntrada dado;
	DadosLinha linhaAux;

	tmpDados = malloc(13 * sizeof(char**));
	// Abre um arquivo TEXTO para LEITURA
	arq = fopen(nomeArquivo, "rt");
	if (arq == NULL) // Se houve erro na abertura
	{
		printf("Problemas na abertura do arquivo\n");
		dado.numReg = 0;
		return dado;
	}
	while (l < 13) {
		// Lê uma linha (inclusive com o '\n')
		result = fgets(Linha, 100, arq); // o 'fgets' lê até 99 caracteres ou até o '\n'
		result = fgets(Linha, 100, arq); // o 'fgets' lê até 99 caracteres ou até o '\n'
		if (result) // Se foi possível ler
		{
			tmpDados[l] = malloc(strlen(Linha) * sizeof(char*));
			strcpy(tmpDados[l], Linha);
			l++;
		}
	}

	dado.numReg = atoi(tmpDados[0]);
	dado.numZon = atoi(tmpDados[1]);

	dado.tamReg = DoubVet1D(dado.numReg);
	linhaAux = LinhaSplit(tmpDados[2], split);
	for (int i = 0; i < dado.numReg; i++) {
		dado.tamReg[i] = atof(linhaAux.linha[i]);
	}

	dado.partReg = IntVet1D(dado.numReg);
	linhaAux = LinhaSplit(tmpDados[3], split);
	for (int i = 0; i < dado.numReg; i++) {
		dado.partReg[i] = atoi(linhaAux.linha[i]);
	}

	dado.fontReg = DoubVet1D(dado.numReg);
	linhaAux = LinhaSplit(tmpDados[4], split);
	for (int i = 0; i < dado.numReg; i++) {
		dado.fontReg[i] = atof(linhaAux.linha[i]);
	}

	dado.zonReg = IntVet1D(dado.numReg);
	linhaAux = LinhaSplit(tmpDados[5], split);
	for (int i = 0; i < dado.numReg; i++) {
		dado.zonReg[i] = atoi(linhaAux.linha[i]);
	}

	dado.stZon = DoubVet1D(dado.numZon);
	linhaAux = LinhaSplit(tmpDados[6], split);
	for (int i = 0; i < dado.numZon; i++) {
		dado.stZon[i] = atof(linhaAux.linha[i]);
	}

	dado.ssZon = DoubVet1D(dado.numZon);
	linhaAux = LinhaSplit(tmpDados[7], split);
	for (int i = 0; i < dado.numZon; i++) {
		dado.ssZon[i] = atof(linhaAux.linha[i]);
	}

	dado.erro = atof(tmpDados[8]);
	dado.ordQuad = atoi(tmpDados[9]);
	dado.cce = atof(tmpDados[10]);
	dado.ccd = atof(tmpDados[11]);
	dado.p = atoi(tmpDados[12]);

	//liberando memória
	for (int i = 0; i < 12; i++)
		free(tmpDados[i]);
	free(tmpDados);

	//fechando o arquivo
	fclose(arq);
	return dado;
}

void Escrever(DadosEntrada dados, DiamodDifference DD) {
	FILE *file;
	file =
			fopen(
					"/home/matheus/Documentos/Projeto_de_Pesquisa/1-Projeto_Oficial/projetoIntegrado/ResultsFiles/Result.txt",
					"w");

	fprintf(file, "//TEMPO-DE-EXECUCAO\n");
	fprintf(file, "%.5f\n", DD.tempoVar);

	fprintf(file, "//FUGA-PARA-DIREITA\n");
	fprintf(file, "%.5f\n", DD.fugaDir);
	fprintf(file, "//FUGA-PARA-ESQUERDA\n");
	fprintf(file, "%.5f\n", DD.fugaEsq);

	fprintf(file, "//FLUXO-ESCALAR\n");
	for (int i = 0; i < DD.dimFi; i += dados.p)
		fprintf(file, "%.2f %.7f\n", DD.x[i], DD.fi[i]);

	fclose(file);
}

DadosLinha LinhaSplit(char str[], char split[]) {
	char *temp = 0;
	char **result = 0;
	unsigned int tamanho = 0;
	DadosLinha tmp;
	temp = strtok(str, split);

	if (temp) {
		result = malloc((tamanho + 1) * sizeof(char**));
		result[tamanho++] = temp;
	}

	while ((temp = strtok(0, split)) != 0) {
		result = realloc(result, (tamanho + 1) * sizeof(char**));
		result[tamanho++] = temp;
	}
	tmp.linha = result;
	tmp.dim = tamanho;
	return tmp;
}

/*
 * Manipulação de vetores e matrizes de números inteiros
 */
int* IntVet1D(int dim1) {
	int *vet1D;

	// alocando o vetor de ponteiros
	vet1D = (int*) malloc(dim1 * sizeof(int));

	return vet1D;
}

int** Intvet2D(int dim1, int dim2) {
	// declaração de variável ponteiro para ponteiro
	int **vet2D;

	// alocando o vetor de ponteiros
	vet2D = (int**) malloc(dim1 * sizeof(int*));

	// alocando cada uma das linhas da matriz
	for (int i = 0; i < dim1; i++)
		vet2D[i] = (int*) malloc(dim2 * sizeof(int));
	return vet2D;
}

void LiberarInt1D(int *vet1D) {
	// liberando o espaço de memória
	// libera o vetor de ponteiros
	free(vet1D);
}

void LiberarInt2D(int **vet2D, int dim1) {
	// liberando o espaço de memória
	// libera cada linha
	for (int i = 0; i < dim1; i++)
		free(vet2D[i]);
	// libera o vetor de ponteiros
	free(vet2D);
}

/*
 * Manipulação de vetores e matrizes de números reais
 */
double* DoubVet1D(int dim1) {
	double *vet1D;

	// alocando o vetor de ponteiros
	vet1D = (double*) malloc(dim1 * sizeof(double));

	return vet1D;
}

double** DoubVet2D(int dim1, int dim2) {
	// declaração de variável ponteiro para ponteiro
	double **vet2D;

	// alocando o vetor de ponteiros
	vet2D = (double**) malloc(dim1 * sizeof(double*));

	// alocando cada uma das linhas da matriz
	for (int i = 0; i < dim1; i++)
		vet2D[i] = (double*) malloc(dim2 * sizeof(double));
	return vet2D;
}

double*** DoubVet3D(int dim1, int dim2, int dim3) {
	// declaração de variável ponteiro para ponteiro
	double ***vet3D;

	// alocando o vetor de ponteiros
	vet3D = (double***) malloc(dim1 * sizeof(double**));

	// alocando cada uma das linhas da matriz
	for (int i = 0; i < dim1; i++) {
		vet3D[i] = (double**) malloc(dim2 * sizeof(double*));
		for (int j = 0; j < dim2; j++)
			vet3D[i][j] = (double*) malloc(dim3 * sizeof(double));
	}
	return vet3D;
}

void LiberarDoub1D(double *vet1D) {
	// liberando o espaço de memória
	// libera o vetor de ponteiros
	free(vet1D);
}

void LiberarDoub2D(double **vet2D, int dim1) {
	// liberando o espaço de memória
	// libera cada linha
	for (int i = 0; i < dim1; i++)
		free(vet2D[i]);
	// libera o vetor de ponteiros
	free(vet2D);
}

void LiberarDoub3D(double ***vet3D, int dim1, int dim2) {
	// liberando o espaço de memória
	// libera cada linha
	for (int i = 0; i < dim1; i++) {
		for (int j = 0; j < dim2; j++)
			free(vet3D[i][j]);
		free(vet3D[i]);
	}
	// libera o vetor de ponteiros
	free(vet3D);
}

/*
 * Gerador de quadraturas de Gauss-Legandre
 */
PLegendre GaussLegendreAbsPes(int n) {
	PLegendre result;
	double pz, p0, p1, dpz, z1;
	double *A = DoubVet1D(n);
	double *P = DoubVet1D(n);

	int Info = 0;
	double Toler = 1e-15;
	int IterMax = 100000;
	int m = (n + 1) / 2;
	double fracn = 1 - (1 - 1 / n) / (8 * n * n);
	double pin = M_PI / (n + 0.5);
	for (int i = 1; i <= m; i++) {
		double z = fracn * cos((i - 0.25) * pin);
		int Iter = 0;
		do {
			Iter++;
			pz = 1;
			p1 = 0;
			for (int j = 0; j < n; j++) {
				p0 = p1;
				p1 = pz;
				pz = ((2 * j + 1) * z * p1 - j * p0) / (j + 1);
			}
			dpz = (n * (p1 - z * pz)) / (1 - z * z);
			z1 = z;
			z = z1 - pz / dpz;
		} while ((Abs(z - z1) >= Toler) && (Iter != IterMax));
		if (Abs(z - z1) <= Toler) {
			A[i - 1] = z;
			A[m + i - 1] = -z;
			P[i - 1] = 2 / ((1 - z * z) * dpz * dpz);
			P[m + i - 1] = P[i - 1];
		} else
			Info++;
	}
	result.mi = A;
	result.w = P;
	//result.Info = Info;
	return result;
}

double Abs(double x) {
	if (x < 0)
		return -x;
	return x;
}

/*
 * Execução dos cálculos do método Diamod Difference DD
 */
DiamodDifference getDD(DadosEntrada e) {
	DiamodDifference dd;
	/*
	 * Preenche os vetores com os dados do domínio
	 */
	int J = 0;
	for (int i = 0; i < e.numReg; i++)
		J += e.partReg[i];
	int N = e.ordQuad;
	dd.dimFi = J + 1;
	//NNR = e.NNR;
	//double Stop = e.erro;
	//TipoContorno = e.CCETipo + e.CCDTipo;
	double *St = DoubVet1D(J);
	double *Ss = DoubVet1D(J);
	double *Q = DoubVet1D(J);
	double *h = DoubVet1D(J);
	double *S = DoubVet1D(J);
	//double *Fi = DoubVet1D(J + 1);
	//double *x = DoubVet1D(J + 1);
	//double *TaxaAbR = DoubVet1D(e.numReg);
	double **psiOld = DoubVet2D(e.ordQuad, J + 1);
	double **psi = DoubVet2D(e.ordQuad, J + 1);
	double tmpErro = 0;
	double erro = 0;
	int k = 0;
	dd.fi = DoubVet1D(J + 1);
	dd.x = DoubVet1D(J + 1);
	dd.x[0] = 0;
	for (int i = 0; i < e.numReg; i++) {
		double temp = e.tamReg[i] / e.partReg[i];
		for (int j = 0; j < e.partReg[i]; j++) {
			h[k] = temp;
			Ss[k] = e.ssZon[e.zonReg[i]];
			St[k] = e.stZon[e.zonReg[i]];
			Q[k] = e.fontReg[i];
			dd.x[k + 1] = dd.x[k] + temp;
			k++;
		}
	}
	/*  Estabelecendo a condição de contorno */
	for (int n = 0; n < e.ordQuad / 2; n++) {
		if (e.cce < 0) {
			psi[n][0] = 0;
			psiOld[n][0] = psi[n][0];
		} else {
			psi[n][0] = e.cce;
			psiOld[n][0] = psi[n][0];
		}
		if (e.ccd < 0) {
			psi[e.ordQuad / 2 + n][J] = 0;
			psiOld[e.ordQuad / 2 + n][0] = psi[e.ordQuad / 2 + n][0];
		} else {
			psi[e.ordQuad / 2 + n][J] = e.ccd;
			psiOld[e.ordQuad / 2 + n][0] = psi[e.ordQuad / 2 + n][0];
		}
	}
	/* Obtendo as quadraturas de Gauss-Legendre */
	PLegendre pl = GaussLegendreAbsPes(e.ordQuad);
	double *Mi = pl.mi;
	double *W = pl.w;
	//#############################################
	double start;
	start = omp_get_wtime();
	dd.numInter = 0;
	do {
		dd.numInter++;
		/* Calculando a fonte de espalhamento */
		for (int j = 0; j < J; j++) {
			S[j] = 0;
			for (int n = 0; n < N; n++)
				S[j] += 0.25 * Ss[j] * (psi[n][j + 1] + psi[n][j]) * W[n];
		}

		/* Varredura para a direita/esquerda */
		for (int m = 0; m < N / 2; m++)
			for (int j = 0; j < J; j++) {
				/* Varredura para a direita */
				psiOld[m][j + 1] = psi[m][j + 1];
				psi[m][j + 1] = ((Mi[m] / h[j] - 0.5 * St[j]) * psi[m][j] + S[j]
						+ Q[j]) / (Mi[m] / h[j] + 0.5 * St[j]);
				/* Varredura para a esquerda */
				psiOld[m + N / 2][J - 1 - j] = psi[m + N / 2][J - 1 - j];
				psi[m + N / 2][J - 1 - j] = ((Mi[m] / h[J - 1 - j]
						- 0.5 * St[J - 1 - j]) * psi[m + N / 2][J - j]
						+ S[J - 1 - j] + Q[J - 1 - j])
						/ (Mi[m] / h[J - 1 - j] + 0.5 * St[J - 1 - j]);
			}

		/*  Estabelecendo a condição de contorno */
		for (int n = 0; n < e.ordQuad / 2; n++) {
			if (e.cce < 0) {
				psiOld[n][0] = psi[n][0];
				psi[n][0] = psi[e.ordQuad / 2 + n][0];
			}
			if (e.ccd < 0) {
				psiOld[e.ordQuad / 2 + n][0] = psi[n][0];
				psi[e.ordQuad / 2 + n][J] = psi[n][J];
			}
		}

		/* Calculando o erro para o critério de parada*/
		tmpErro = 0;
		erro = 0;
		for (int m = 0; m < N; m++)
			for (int j = 0; j < J; j++) {
				tmpErro = Abs(psi[m][j] - psiOld[m][j]) / psi[m][j];
				if (tmpErro > erro)
					erro = tmpErro;
			}
	} while (erro > e.erro);
	dd.tempoVar = omp_get_wtime() - start;
	/* Calcula o fluxo escalar */
	for (int j = 0; j <= J; j++) {
		dd.fi[j] = 0;
		for (int m = 0; m < N; m++)
			dd.fi[j] += psi[m][j] * W[m];
	}

	dd.fugaDir = 0;
	dd.fugaEsq = 0;

	for (int m = 0; m < e.ordQuad / 2; m++) {
		dd.fugaDir += Mi[m] * psi[m][J] * W[m];
	}

	for (int m = e.ordQuad / 2; m < e.ordQuad; m++) {
		dd.fugaEsq += Mi[m] * psi[m][0] * W[m];
	}

	/*
	 double taxaABS[] = new double[numRegioes];
	 double sigmaAbs;

	 for (int i = 0; i < e.numReg; i++) {
	 sigmaAbs = St[i] - Ss[i];

	 for (int j = 0; j < J; j++) {
	 taxaABS[i] += sigmaAbs * h[i] * fluxoEscalarMedio()[j];
	 }
	 }
	*/

	/*
	int i = 0;
	double TaxaAbR[], TaxaAbT[];

	for (int r = 0; r < e.numReg; r++) {
		for (int j = i; j < i + NNR[r]; j++) {
			TaxaAbR[r] += 0.5 * (dd.fi[j] + dd.fi[j + 1]);
		}
		TaxaAbR[r] *= (e.stZon[i] - e.ssZon[i]) * h[i];
		i += NNR[r];
	}
	TaxaAbT = TaxaAbR.Sum;
	*/

	return dd;
}

DiamodDifference getDDParallel(DadosEntrada e) {
	DiamodDifference dd;
	/*
	 * Preenche os vetores com os dados do domínio
	 */
	int J = 0;
	for (int i = 0; i < e.numReg; i++)
		J += e.partReg[i];
	int N = e.ordQuad;
	dd.dimFi = J + 1;
	//NNR = e.NNR;
	//double Stop = e.erro;
	//TipoContorno = e.CCETipo + e.CCDTipo;
	double *St = DoubVet1D(J);
	double *Ss = DoubVet1D(J);
	double *Q = DoubVet1D(J);
	double *h = DoubVet1D(J);
	double *S = DoubVet1D(J);
	//double *Fi = DoubVet1D(J + 1);
	//double *x = DoubVet1D(J + 1);
	//double *TaxaAbR = DoubVet1D(e.numReg);
	double **psiOld = DoubVet2D(e.ordQuad, J + 1);
	double **psi = DoubVet2D(e.ordQuad, J + 1);
	double tmpErro = 0;
	double erro = 0;
	int k = 0;
	dd.fi = DoubVet1D(J + 1);
	dd.x = DoubVet1D(J + 1);
	dd.x[0] = 0;
	for (int i = 0; i < e.numReg; i++) {
		double temp = e.tamReg[i] / e.partReg[i];
		for (int j = 0; j < e.partReg[i]; j++) {
			h[k] = temp;
			Ss[k] = e.ssZon[e.zonReg[i]];
			St[k] = e.stZon[e.zonReg[i]];
			Q[k] = e.fontReg[i];
			dd.x[k + 1] = dd.x[k] + temp;
			k++;
		}
	}
	/*  Estabelecendo a condição de contorno */
	for (int n = 0; n < e.ordQuad / 2; n++) {
		if (e.cce < 0) {
			psi[n][0] = 0;
			psiOld[n][0] = psi[n][0];
		} else {
			psi[n][0] = e.cce;
			psiOld[n][0] = psi[n][0];
		}
		if (e.ccd < 0) {
			psi[e.ordQuad / 2 + n][J] = 0;
			psiOld[e.ordQuad / 2 + n][0] = psi[e.ordQuad / 2 + n][0];
		} else {
			psi[e.ordQuad / 2 + n][J] = e.ccd;
			psiOld[e.ordQuad / 2 + n][0] = psi[e.ordQuad / 2 + n][0];
		}
	}
	/* Obtendo as quadraturas de Gauss-Legendre */
	PLegendre pl = GaussLegendreAbsPes(e.ordQuad);
	double *Mi = pl.mi;
	double *W = pl.w;
	//#############################################
	double start;
	omp_set_num_threads(MAX_THREADS);
	start = omp_get_wtime();
	dd.numInter = 0;
	do {
		dd.numInter++;
		/* Calculando a fonte de espalhamento */
#pragma omp parallel for
		for (int j = 0; j < J; j++) {
			S[j] = 0;
			for (int n = 0; n < N; n++)
				S[j] += 0.25 * Ss[j] * (psi[n][j + 1] + psi[n][j]) * W[n];
		}

		/* Varredura para a direita/esquerda */
#pragma omp sections
		{
#pragma omp section
			{
#pragma omp parallel for
				for (int m = 0; m < N / 2; m++)
					for (int j = 0; j < J; j++) {
						/* Varredura para a direita */
						psiOld[m][j + 1] = psi[m][j + 1];
						psi[m][j + 1] = ((Mi[m] / h[j] - 0.5 * St[j])
								* psi[m][j] + S[j] + Q[j])
								/ (Mi[m] / h[j] + 0.5 * St[j]);
					}
			}
#pragma omp section
			{
#pragma omp parallel for
				for (int m = 0; m < N / 2; m++)
					for (int j = 0; j < J; j++) {
						/* Varredura para a esquerda */
						psiOld[m + N / 2][J - 1 - j] =
								psi[m + N / 2][J - 1 - j];
						psi[m + N / 2][J - 1 - j] = ((Mi[m] / h[J - 1 - j]
								- 0.5 * St[J - 1 - j]) * psi[m + N / 2][J - j]
								+ S[J - 1 - j] + Q[J - 1 - j])
								/ (Mi[m] / h[J - 1 - j] + 0.5 * St[J - 1 - j]);
					}
			}
		} // Aqui há uma barreira implícita

		/*  Estabelecendo a condição de contorno */
#pragma omp parallel for
		for (int n = 0; n < e.ordQuad / 2; n++) {
			if (e.cce < 0) {
				psiOld[n][0] = psi[n][0];
				psi[n][0] = psi[e.ordQuad / 2 + n][0];
			}
			if (e.ccd < 0) {
				psiOld[e.ordQuad / 2 + n][0] = psi[n][0];
				psi[e.ordQuad / 2 + n][J] = psi[n][J];
			}
		}

		/* Calculando o erro para o critério de parada*/
		tmpErro = 0;
		erro = 0;
		for (int m = 0; m < N; m++)
			for (int j = 0; j < J; j++) {
				tmpErro = Abs(psi[m][j] - psiOld[m][j]) / psi[m][j];
				if (tmpErro > erro)
					erro = tmpErro;
			}
	} while (erro > e.erro);
	dd.tempoVar = omp_get_wtime() - start;
	/* Calcula o fluxo escalar */
	for (int j = 0; j <= J; j++) {
		dd.fi[j] = 0;
		for (int m = 0; m < N; m++)
			dd.fi[j] += psi[m][j] * W[m];
	}

	dd.fugaDir = 0;
	dd.fugaEsq = 0;

#pragma omp parallel
	{
#pragma omp sections
		{
#pragma omp section
			{
#pragma omp parallel for
				for (int m = 0; m < e.ordQuad / 2; m++) {
					dd.fugaDir += Mi[m] * psi[m][J] * W[m];
				}
			}

#pragma omp section
			{
#pragma omp parallel for
				for (int m = e.ordQuad / 2; m < e.ordQuad; m++) {
					dd.fugaEsq += Mi[m] * psi[m][0] * W[m];
				}
			}
		}
	}

	return dd;
}
