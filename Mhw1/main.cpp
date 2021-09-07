#include <stdexcept>
#include "readfile.h"
#include "eval.h"

using namespace std;

int main()
{
  vector<Atom> Atoms;
  string fname = "example.txt";
  try
  {
    Readatomsformfile(Atoms, fname);
  }
  catch (invalid_argument &e)
  {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }

  // cout << E_LJ(Atoms)<< endl;
  // arma::mat Force(3, 3);
  // F_LJ_fd(Force, Atoms, 0.0001);
  // Force.print("Force");
  // F_LJ_fd(Force, Atoms, 0.00001);
  // Force.print("Force");

  vector<Atom> opt_Atoms = Atoms;
  // Steepest_descend( opt_Atoms, Atoms, fdstepsize, double thresh)
  Steepest_descend( opt_Atoms, Atoms,  0.0001, 1e-2);
  for (auto x : opt_Atoms)
    cout << x << std::endl;

  return EXIT_SUCCESS;
}
