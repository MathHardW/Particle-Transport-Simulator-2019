//===========================STRUCT'S===================================
typedef struct {
	double *AutoValor;
	double **AutoVetor;
} AutosVV;
//======================================================================

//===========================CABEÇALHOS=================================

double *multiplicarMatrizes(int ORDER, double *A, double *B, double *C);
double Abs(double x);
double *inverse(double *A, int N);

AutosVV *SetAutosVV(int OQ, int ZM);
AutosVV *CalcAutosVV(double ***matrizMR, int OQ, int ZM);
double *print_eigenvalues( char* desc, int n, double* wr, double* wi);
double **print_eigenvectors( char* desc, int n, double* wi, double* v,
							 int ldv);

int     *AlocarInteger1D(int Dim1);
double  *AlocarDouble1D (int Dim1);
double **AlocarDouble2D (int Dim1, int Dim2);
double ***AlocarDouble3D(int Dim1, int Dim2, int Dim3);
//======================================================================

//=============================ROTINAS==================================

// LU decomoposition of a general matrix
extern void dgetrf_(int* M, int *N, double* A, int* lda,
					int* IPIV, int* INFO);
// generate inverse of a matrix given its LU decomposition
extern void dgetri_(int* N, double* A, int* lda, int* IPIV,
					double* WORK, int* lwork, int* INFO);
//DGEEV computes the eigenvalues and eigenvectors for GE matrices
extern void dgeev_( char* jobvl, char* jobvr, int* n, double* a,
					int* lda, double* wr, double* wi, double* vl,
					int* ldvl, double* vr, int* ldvr, double* work,
					int* lwork, int* info );
//======================================================================

//=============================FUNÇÕES==================================

//GAUSS LEGENDRE========================================================
typedef struct{
    double *mi;   // define o campo direção mi
    double *w;   // define o peso w
} PLegendre;    // Define o nome do novo tipo criado

PLegendre GaussLegendreAbsPes(int n);

PLegendre GaussLegendreAbsPes(int n){
	PLegendre result;
	double pz, p0, p1, dpz, z1;
	double *A = AlocarDouble1D(n);
	double *P = AlocarDouble1D(n);

	int Info = 0;
	double Toler = 1e-15;
	int IterMax = 100000;
	int m = (n + 1) / 2;
	double fracn = 1 - (1 - 1 / n) / (8 * n * n);
	//double M_PI = 3.14159265358979323846264338327;
	double pin = M_PI / (n + 0.5);
	for (int i = 1; i <= m; i++)
	{
		double z = fracn * cos((i - 0.25) * pin);
		int Iter = 0;
		do
		{
			Iter++;
			pz = 1;
			p1 = 0;
			for (int j = 0; j < n; j++)
			{
				p0 = p1;
				p1 = pz;
				pz = ((2 * j + 1) * z * p1 - j * p0) / (j + 1);
			}
			dpz = (n * (p1 - z * pz)) / (1 - z * z);
			z1 = z;
			z = z1 - pz / dpz;
		} while ((Abs(z - z1) >= Toler) && (Iter != IterMax));
		if (Abs(z - z1) <= Toler)
		{
			A[i - 1] = z;
			A[m + i - 1] = -z;
			P[i - 1] = 2 / ((1 - z * z) * dpz * dpz);
			P[m + i - 1] = P[i - 1];
		}
		else
			Info++;
	}
	result.mi = A;
	result.w = P;
	//result.Info = Info;
	return result;
}
//======================================================================

//MULTIPLICAR MATRIZES==================================================
double *multiplicarMatrizes(int ORDER, double *A, double *B, double *C) {
	//DIMENSÕES .: A[N][P], B[P][M], C[N][M];
    int Ndim, Pdim, Mdim;
    Ndim = Pdim = Mdim = ORDER;

    /* multiply the matrices */	double tmp;
    for (int i = 0; i < Ndim; i++){
        for (int j = 0; j < Mdim; j++){
            tmp = 0.0;
            for(int k = 0; k < Pdim; k++){
                /* C(i,j) = sum(A(i,k) * B(k,j)) */
                tmp += (*(A+(i*Ndim+k))) * (*(B+(k*Pdim+j)));
            }
            *(C+(i*Ndim+j)) = tmp;
        }
    }
}
//======================================================================

//CALCULAR INVERSA DA MATRIZ============================================
double *inverse(double *A, int N){
    int *IPIV = (int*)malloc( (N+1)*sizeof(int) );
    int LWORK = N*N;
    double *WORK = (double*)malloc( LWORK*sizeof(double) );
    int INFO;

    dgetrf_(&N,&N,A,&N,IPIV,&INFO);
    dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

    free(IPIV);
    free(WORK);
}
//======================================================================

//AUTOVALOR=============================================================
AutosVV *SetAutosVV(int OQ, int ZM){
	AutosVV *autosVV = (AutosVV*)calloc(ZM*OQ*OQ*OQ, sizeof(AutosVV));

		for(int z=0; z<ZM; z++){
			autosVV[z].AutoValor = AlocarDouble1D(OQ);
			autosVV[z].AutoVetor = AlocarDouble2D(OQ, OQ);
		}

	return autosVV;
}

AutosVV *CalcAutosVV(double ***matriz, int OQ, int ZM) {
        AutosVV *autosVV = SetAutosVV(OQ, ZM);

        //VARIÁVEIS INTEGERS============================================
        int n = OQ, lda = n, ldvl = n, ldvr = n, info, lwork, cont=0;
        //VARIÁVEIS DOUBLES=============================================
        double wkopt;
        double* work;
        //VARIÁVEIS VETORES=============================================
        double wr[n], wi[n], vl[ldvl*n], vr[ldvl*n];
        double* a = AlocarDouble1D(ZM*OQ*OQ);


		for(int z=0; z<ZM; z++){
			for(int l=0; l<OQ; l++){
				for(int c=0; c<OQ; c++){
					a[cont] = matriz[z][l][c];
					cont++;
				}
			}

				/*EXECUÇÕES DA CLASSE ORIGINAL*/{
					//printf( " DGEEV Example Program Results\n" );
					/* Query and allocate the optimal workspace */
					lwork = -1;
					dgeev_( "Vectors", "Vectors", &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
					 &wkopt, &lwork, &info );
					lwork = (int)wkopt;
					work = (double*)malloc( lwork*sizeof(double) );
					/* Solve eigenproblem */
					dgeev_( "Vectors", "Vectors", &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
					 work, &lwork, &info );
					/* Check for convergence */
					if( info > 0 ) {
							printf( "The algorithm failed to compute eigenvalues.\n" );
							exit( 1 );
					}
				}

			//CAPTURANDO E IMPRIMINDO AUTOVALOR===================================
			autosVV[z].AutoValor = print_eigenvalues( "Eigenvalues", n, wr, wi );
			//====================================================================

			//CAPTURANDO E IMPRIMINDO AUTOVETOR DA ESQUERDA======================
			autosVV[z].AutoVetor = print_eigenvectors( "Left eigenvectors", n, wi, vl, ldvl);
			//====================================================================

			cont = 0;
		}


        /* LIBERANDO MEMÓRIA */
        free( (void*)work );

        return autosVV;
}

//======================================================================
/* Auxiliary routine: printing eigenvalues */
double *print_eigenvalues( char* desc, int n, double* wr, double* wi) {
	double *AutoValor = AlocarDouble1D(n);

    int j;
    //printf( "\n %s\n", desc );

	for( j = 0; j < n; j++ ) {
		if( wi[j] == (double)0.0 ) {
			AutoValor[j] = wr[j];
			//printf( " %f", wr[j]);
		} else {
			AutoValor[j] = wr[j];
			//printf( " (%f,%f)", wr[j], wi[j] );
		}
	}
	//printf( "\n" );

	return AutoValor;
}

/* Auxiliary routine: printing eigenvectors */
double **print_eigenvectors( char* desc, int n, double* wi, double* v, int ldv) {
   double  *AutoVetorTemp = AlocarDouble1D(n*n);
   double **AutoVetor = AlocarDouble2D(n,n);

   int i, j, cont=0;
   //printf( "\n %s\n", desc );
   for( i = 0; i < n; i++ ) {
      j = 0;
      while( j < n ) {
         if( wi[j] == (double)0.0 ) {
			AutoVetorTemp[cont] = v[i+j*ldv];
            //printf( " %f", v[i+j*ldv]);
            j++; cont++;
         } else {
			AutoVetorTemp[cont] = v[i+j*ldv];
			AutoVetorTemp[cont+1] = v[i+j*ldv];
            //printf( " (%f,%f)", v[i+j*ldv], v[i+(j+1)*ldv] );
            //printf( " (%f,%f)", v[i+j*ldv], -v[i+(j+1)*ldv] );
            j += 2; cont +=2;
         }
      }
      //printf( "\n" );
   }

	cont=0;
	for(int j=0; j<n; j++){
		for(int x=0; x<n; x++){
			AutoVetor[j][x] = AutoVetorTemp[cont];
			cont++;
		}
	}

	free(AutoVetorTemp);
	return AutoVetor;
}
//======================================================================

//MÉTODO MÓDULO=========================================================
double Abs(double x){
	if(x < 0)
		return -x;
	return x;
}
//======================================================================
