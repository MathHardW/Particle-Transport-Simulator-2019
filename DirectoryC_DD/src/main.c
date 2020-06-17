#include "DDHeader.h"

int main(int argc, char *argv[]) {
	printf("start\n\n");

	int metodo = atoi(argv[2]);

	DadosEntrada dados = Ler(argv[1]);

	if (metodo == 0) {
		DiamodDifference DD = getDD(dados);
		//for (int i = 0; i < DD.dimFi; i += dados.p) printf("x = %.2f\tFluxo Escalar serial= %.7f\n", DD.x[i], DD.fi[i]);

		printf(
				"Diamond Diference Serial---------------------------------\nIterações %i, Tempo de varredura igual a %.5f segundos.\n\n",
				DD.numInter, DD.tempoVar);

		Escrever(dados, DD);
	} else {
		DiamodDifference DDParallel = getDDParallel(dados);
		//for (int i = 0; i < DDParallel.dimFi; i += dados.p) printf("x = %.2f\tFluxo Escalar paralelo= %.7f\n", DDParallel.x[i], DDParallel.fi[i]);

		printf(
				"Diamond Diference Paralelo---------------------------------\nIterações %i, Tempo de varredura igual a %.5f segundos.\n\n",
				DDParallel.numInter, DDParallel.tempoVar);

		Escrever(dados, DDParallel);
	}

	printf("finish");

	return 0;
}
