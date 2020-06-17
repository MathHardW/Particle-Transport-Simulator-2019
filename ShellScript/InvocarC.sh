#!/bin/bash

cd /home/matheus/Documentos/Projeto_de_Pesquisa/1-Projeto_Oficial/projetoIntegrado/DirectoryC_DD/src 
gcc -o main main.c DDHeader.h DD.c -I lib -lopenblas -fopenmp -lm 
./main /home/matheus/Documentos/Projeto_de_Pesquisa/1-Projeto_Oficial/projetoIntegrado/DataFiles/DOMINGUES 1