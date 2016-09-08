//Solver headers
//solver type; if not define SOLVER_MT, using superLU 
//#define SOLVER_MT
#ifdef SOLVER_MT
#include "../solver/superLUMT_solver.h"
#else
#include "../solver/superLU_solver.h"
#endif
