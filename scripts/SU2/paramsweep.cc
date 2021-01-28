#include "su2.h"

int main(int argc, char* argv[]) {
  if(argc < 2) 
    { 
    printfln("Usage: %s inputfile",argv[0]); 
    return 0; 
    }
  auto input = InputGroup(argv[1],"input");
  auto Nmin = input.getInt("Nmin", 10);
  auto Nmax = input.getInt("Nmax", 50);
  auto dN   = input.getInt("dN", 10);
  auto amin = input.getReal("amin", 0.2);
  auto amax = input.getReal("amax", 1);
  auto da   = input.getReal("da", 0.2);
  auto mmin = input.getReal("mmin", 0);
  auto mmax = input.getReal("mmax", 1);
  auto dm   = input.getReal("dm", 0.25);
  auto gmin = input.getReal("gmin", 0);
  auto gmax = input.getReal("gmax", 1);
  auto dg   = input.getReal("dg", 0.25);
  auto nstates   = input.getInt("nstates", 1);
  auto weight   = input.getReal("weight", 100);

  int N;
  double a;
  double m;
  double g;
  auto fpath = "./run";
  FILE *runfile;
  for (N = Nmin; N <= Nmax; N+=dN) {
    for (a = amin; a <= amax; a+=da) {
      for (m = mmin; m <= mmax; m+=dm) {
	for (g = gmin; g <= gmax; g+=dg) {
	  runfile = std::fopen(fpath, "w");
	  std::fprintf(runfile, "input\n{\nN=%d\na=%.12f\nm=%.12f\ng=%.12f\nnstates=%d\nweight=%.12f\n}", N, a, m, g, nstates, weight);
	  std::fclose(runfile);
      
      su2(fpath);
      // remove file
      std::remove(fpath);      
	}
      }
    }
  }
  return 0;
}
