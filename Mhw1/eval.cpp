#include "eval.h"
#include <math.h>
#include <cassert>
#include <armadillo>

using namespace std;
double e_Au = 2.951;
double r_Au =5.29; 
// double r_Au =2.0; 
// double e_Au = 2;

void ConvertAtomstoCoor(const vector<Atom> &Atoms, arma::mat &Coor){
    for(size_t k = 0; k < Atoms.size(); k++)
        for (size_t j = 0; j < 3; j++)
            Coor(j,k) = Atoms[k].r[j];
}

void ConvertCoortoAtoms(vector<Atom> &Atoms, const arma::mat &Coor){
    for(size_t k = 0; k < Atoms.size(); k++)
        for (size_t j = 0; j < 3; j++)
            Atoms[k].r[j] = Coor(j,k);
}

double E_LJ(const std::vector<Atom> &Atoms)
{
    double E = 0.0;
    for (size_t k = 0; k < Atoms.size(); k++){
        Atom atom_k = Atoms[k];
        for (size_t j = k + 1; j < Atoms.size(); j++)
        {
            Atom atom_j = Atoms[j];
            double R2 = (atom_j.r[0] - atom_k.r[0]) * (atom_j.r[0] - atom_k.r[0])
                 + (atom_j.r[1] - atom_k.r[1]) * (atom_j.r[1] - atom_k.r[1]) 
                 + (atom_j.r[2] - atom_k.r[2]) * (atom_j.r[2] - atom_k.r[2]);
            double r2overR2 =  r_Au * r_Au / R2;
            double E_temp = pow(r2overR2, 6) - 2.0 * pow(r2overR2, 3);
            E += E_temp * e_Au;
            // cout << E_temp * e_Au << endl;
        }
    }
    return E;
}

void F_LJ_fd(arma::mat &force, const vector<Atom> &Atoms, double stepsize)
{
    vector<Atom> Atoms_forward = Atoms;
    vector<Atom> Atoms_backward = Atoms;
    assert(force.n_rows == 3 && force.n_cols == Atoms.size());
    for(size_t k = 0; k < Atoms.size(); k++)
        for (size_t j = 0; j < 3; j++){
            Atoms_forward = Atoms;
            Atoms_backward = Atoms;
            Atoms_forward[k].r[j] = Atoms_forward[k].r[j] + stepsize;
            Atoms_backward[k].r[j] = Atoms_backward[k].r[j] - stepsize;
            // cout << Atoms_forward[k].r[j] << endl;
            // cout << Atoms_backward[k].r[j] << endl;
            double E_forward = E_LJ(Atoms_forward);
            double E_backward = E_LJ(Atoms_backward);
            // cout << E_forward << endl;
            // cout << E_backward << endl;
            force(j,k) = -(E_forward - E_backward) /2.0 /stepsize;
        }
}

void Steepest_descend( vector<Atom> &opt_Atoms, const vector<Atom> &Atoms, double fdstepsize, double thresh){

    double search_stepsize = 0.3;
    double E_old = E_LJ(Atoms), E_new;
    vector<Atom> new_Atoms = Atoms;
    arma::mat new_point(3, Atoms.size());
    arma::mat old_point(3, Atoms.size());
    arma::mat Force(3, Atoms.size());
    ConvertAtomstoCoor(Atoms, old_point);

    // evaluate the Force at starting point
    F_LJ_fd(Force, Atoms, fdstepsize);
    Force.print("Force");
    cout << "Force norm: "<< arma::norm(Force, "fro") << endl;
    cout << "Energy: "<< E_LJ(Atoms) << endl;
    int count = 0;
    cout << endl << "Start opt" << endl << endl;
    while (arma::norm(Force, "fro") > thresh && count < 1e2){
        // calculate new point position
        new_point = old_point + search_stepsize * Force / arma::norm(Force, 2);
        ConvertCoortoAtoms(new_Atoms, new_point);
        E_new = E_LJ(new_Atoms);
        cout << "count "<< count  << " Energy: "<< E_new << endl;
        // the step makes function evaluation lower - it is a good step. what do you do?
        if (E_new < E_old){
            old_point = new_point;
            E_old = E_new;
            // refresh Force
            F_LJ_fd(Force, new_Atoms, fdstepsize);
            search_stepsize *= 1.2;
            new_point.print("new_point");
            cout << "Force norm: "<< arma::norm(Force, "fro") << endl;
        }
        else // the step makes function evaluation higher - it is a bad step. what do you do?
            search_stepsize /= 2;
        count += 1;

    }

    cout << "Total iterations: " << count << endl;
    ConvertCoortoAtoms(opt_Atoms, old_point);
}
