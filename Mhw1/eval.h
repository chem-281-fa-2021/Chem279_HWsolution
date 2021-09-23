#ifndef EVAL_H
#define EVAL_H

#include "readfile.h"
#include <armadillo>

double E_LJ(const std::vector<Atom> &Atoms);

void F_LJ_fd(arma::mat &force, const std::vector<Atom> &Atoms, double stepsize);

void Steepest_descend( std::vector<Atom> &opt_Atoms, const std::vector<Atom> &Atoms, double fdstepsize, double thresh);

#endif // EVAL_H