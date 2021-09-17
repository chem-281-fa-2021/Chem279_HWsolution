#include <stdlib.h>
#include <stdexcept>
#include <armadillo>
#include <stdio.h>
#include <Shell.h>
#include <Numerical.h>

using namespace std;

int main(int argc, char* argv[])
{
  if (argc !=2)
  {
  printf("usage hw1 filename, for example hw1 example.txt");
  return EXIT_FAILURE;
  }
  string fname(argv[1]);
  Shell sh1,sh2;
  try
  {
    ReadShellparameter(sh1, sh2, fname);
  }
  catch (invalid_argument &e)
  {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }

  // sh1.printinfo();
  // sh2.printinfo();

  // printf("This shell has %d functions.\n", sh2.dim_l());
  // const double s1s1 = Overlap_onedim(0.0, 0.0, 1.0, 1.0, 0, 0);
  // printf("Analytical Integral %1.17e for s1s1\n", s1s1);
  // const double s1p1 = Overlap_onedim(0.0, 0.0, 1.0, 1.0, 0, 1);
  // printf("Analytical Integral %1.17e for s1p1\n", s1p1);
  // const double s1s2 = Overlap_onedim(1.0, 0.0, 1.0, 1.0, 0, 0);
  // printf("Analytical Integral %1.17e for s1s2\n", s1s2);
  // const double s1p2 = Overlap_onedim(0.0, 1.0, 1.0, 1.0, 0, 1);
  // printf("Analytical Integral %1.17e for s1p2\n", s1p2);

  int dim1 = sh1.dim_func(), dim2 = sh2.dim_func();
  arma::mat Overlap_matrix(dim1, dim2,arma::fill::zeros);
  Eval_Ov(Overlap_matrix, sh1, sh2);
  Overlap_matrix.print();

  return EXIT_SUCCESS;
}
