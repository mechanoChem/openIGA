/**
* superLU_MT solver declartion
*/

#ifndef SUPERLUMT_H_
#define SUPERLUMT_H_
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pdsp_defs.h>

int luSolver(double *aTemp, int *asubTemp, int *xaTemp, double *rhsTemp, int m, int n, int nnz, double * _sol, int _solverInfo);

#endif