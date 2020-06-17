//===========================CABEÇALHOS=================================

// 1 DIMENSÃO===========================================================
int *AlocarInteger1D(int Dim1);
double *AlocarDouble1D(int Dim1);
// 2 DIMENSÕES==========================================================
double **AlocarDouble2D(int Dim1, int Dim2);
// 3 DIMENSÕES==========================================================
double ***AlocarDouble3D(int Dim1, int Dim2, int Dim3);
//======================================================================



//=============================FUNÇÕES==================================

// 1 DIMENSÃO===========================================================
int *AlocarInteger1D(int Dim1){
	int *vet1D;

	// alocando o vetor de ponteiros
	vet1D = (int*)malloc(Dim1 * sizeof(int));

	return vet1D;
}

double *AlocarDouble1D(int Dim1){
	double *vet1D;

	// alocando o vetor de ponteiros
	vet1D = (double*)calloc(Dim1, sizeof(double));

	return vet1D;
}

// 2 DIMENSÕES==========================================================
double **AlocarDouble2D(int Dim1, int Dim2){
	double **matriz = malloc(sizeof(double*) * Dim1);

	for(int i=0; i<=Dim2; i++){
		matriz[i] = malloc(sizeof (double) * Dim2);
	}

	return matriz;
}

// 3 DIMENSÕES==========================================================
double ***AlocarDouble3D(int Dim1, int Dim2, int Dim3){
	double ***matriz = malloc(sizeof(double**) * Dim1);

	for(int i=0; i<=Dim2; i++){
		matriz[i] = malloc(sizeof (double*) * Dim2);
		for(int j=0; j<=Dim3; j++){
			matriz[i][j] = malloc(sizeof (double) * Dim3);
		}
	}

	return matriz;
}
//======================================================================
