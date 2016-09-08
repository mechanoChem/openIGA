/*
* superLU solver declartion
*/

#ifndef SUPERLU_H_
#define SUPERLU_H_
#include "slu_ddefs.h"

int luSolver(double *aTemp, int *asubTemp, int *xaTemp, double *rhsTemp, int m, int n, int nnz, double * _sol, int _solverInfo);

#endif