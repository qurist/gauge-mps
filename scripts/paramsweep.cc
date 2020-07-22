#include "schwinger.h"

int main(int argc, char* argv[]) {
  if(argc < 2) 
    { 
    printfln("Usage: %s inputfile",argv[0]); 
    return 0; 
    }
  auto input = InputGroup(argv[1],"input");
  auto Nmin = input.getInt("Nmin", 20);
  auto Nmax = input.getInt("Nmax", 300);
  auto dN   = input.getInt("dN", 20);
  auto xmin = input.getReal("xmin", 50);
  auto xmax = input.getReal("xmax", 1000);
  auto dx   = input.getReal("dx", 50);
  auto NE = input.getInt("NE", 5);
  auto mg = input.getReal("mg",0.25);
  auto gauss = input.getInt("gauss", 0);
  auto Q     = input.getInt("Q",0);
  auto q     = input.getInt("q",1);

  int N;
  double x;
  auto fpath = "./run";
  FILE *runfile;
  for (N = Nmin; N <= Nmax; N+=dN) {
    for (x = xmin; x <= xmax; x+=dx) {
      runfile = std::fopen(fpath, "w");
      std::fprintf(runfile, "input\n{\nN=%d\nx=%.12f\nNE=%d\nmg=%.12f\ngauss=%i\nQ=%d\nq=%d\n}", N, x, NE, mg, gauss, Q, q);
      std::fclose(runfile);
      
      schwinger(fpath);
      // remove file
      std::remove(fpath);      
    }
  }
  // for x in $(seq 10 10 200); do for N in $(seq 10 10 300); do printf "input\n{\nN=$N\nNE=5\nx=$x\nmu=5\nDmax=200\n}" >> ../data/inputs/input_N$((N))_x$((x))_Emax5_mu5_Dmax200; done; done
  return 0;
}
