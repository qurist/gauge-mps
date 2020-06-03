#include "itensor/all.h"
#include "ladder.h"

using namespace itensor;

using Schwinger = MixedSiteSet<LadderSite, SpinHalfSite>;


// std::tuple<Real,MPS> dodmrg(Args const& args = Args::global()){
//   N = args.getInt("N", 100);
//   NE = args.getInt("NE", 5);
//   x  = args.getReal("x", 100);
//   mu = args.getReal("mu", 5);
//   J  = args.getReal("J", 1);
//   Dmax  = args.getInt("Dmax", 200);

//   return 0;
// }


//int main(int argc, char* argv[]){
int schwinger(const char *inputfile, Args const& args = Args::global()){
  // if(argc < 2) 
  //   { 
  //   printfln("Usage: %s inputfile",argv[0]); 
  //   return 0; 
  //   }
  //auto inputfile = args.getString("input=");
  auto input = InputGroup(inputfile,"input");
  auto N = input.getInt("N", 30);
  auto NE = input.getInt("NE", 5);
  auto x  = input.getReal("x", 100);
  auto mg = input.getReal("mg", 0.25);
  auto mu = 2*sqrt(x)*mg;
  auto gauss = input.getInt("gauss",0);

  auto sites = Schwinger(2*N+1, {"ConserveSz=",true, "ConserveNb=", false, "MaxOcc=",2*NE});//SiteSet(N,4*NE+2);
  auto ampo = AutoMPO(sites);

  auto J = 1;
  // auto mu = 5*J;
  //  auto x = 200*J;
  auto l = 0*J;
  auto s = "";
  auto sign = 1;
  
  for(int j=2; j <= 2*N; j+=2)
    {
      // Fix site-dependent variables
      if (j%4==2) {
	s = "projUp";
	sign = -1;
      }
      else if (j%4==0) {
	s = "projDn";
	sign = 1;
      }
      // Pure Gauge Term
      ampo += J, "N",j+1,"N",j+1;
      ampo += -2*NE*J, "N",  j+1;
      ampo += NE*NE*J, "Id", j+1;
      // Matter Term
      ampo += mu,s,j;
      // Matter-Gauge Term
      if (j!=2*N) {
	ampo += x,"S+",j,"N+",j+1,"S-",j+2;
	ampo += x,"S+",j+2,"N-",j+1,"S-",j;
      }

      //Gauge Lagrangian sum_i lambda * [ Gi^2 ]
      ampo += l,"N",j-1, "N",j-1;
      ampo += l,"N",j+1, "N",j+1;
      ampo += -2*l,"N",j-1, "N",j+1;      
      ampo += l, s,j;
      ampo += 2*sign*l, s, j, "N",j+1;
      ampo += -2*sign*l, s, j, "N", j-1;
    }
  
  auto H = toMPO(ampo);

  auto sweeps     = Sweeps(12);
  sweeps.noise()  = 1E-6,1E-8,1E-10,1E-10,1E-10,1E-12,1E-12,1E-12,1E-12,0,0,0;
  sweeps.maxdim() = 10,20,40,40,80,80,160,160,320,320,640,640;
  sweeps.cutoff() = 1E-10;
  double dE       = 1E-6;


  auto state = InitState(sites);
  for(auto j : range1(2*N+1))
    {
      if(j%4==2) state.set(j,"Dn");
      else if(j%4==0) state.set(j,"Up");
      else state.set(j,str(NE));
    }
  auto psiInit0 = MPS(state);

  auto [energy,psi0] = dmrg(H, psiInit0, sweeps, {"Quiet=",true, "EnergyErrgoal=",1E-6});
  printfln("E(N,x) N x NE mu = %.12f %d %.12f %d %.12f\n",energy/(2*N*x), N, x, NE, mu);

  // // Schmidt spectrum at half-cut
  // auto b = N/2;
  // psi0.position(b); 
  // auto r = leftLinkIndex(psi0,b);
  // auto t = siteIndex(psi0,b);
  // auto [U,S,V] = svd(psi0(b),{r,t});
  // auto u = commonIndex(U,S);
  // for(auto n : range1(dim(u)))
  //   {
  //     printfln("%.12f", elt(S,n,n));
  //   }
  
  // Gauss Law
  if (gauss) {
    for(int j=2; j <= 2*N; j+=2)
      {
        // Site-dependent variables
        if (j%4==2) {
  	s = "projUp";
  	sign = -1;
        }
        else if (j%4==0) {
  	s = "projDn";
  	sign = 1;
        }
        // Construct Gauge operator G_i
        ampo = AutoMPO(sites);
        ampo += -1,"N",j-1;
        ampo += 1,"N",j+1;
        //ampo += -2*NE,"Id",j+1;
        ampo += sign, s,j; 
        auto Gi = toMPO(ampo);
        
        // Construct Gauge-squared G_i^2
        ampo = AutoMPO(sites); 
        ampo += 1,"N",j-1, "N",j-1;
        ampo += 1,"N",j+1, "N",j+1;
        ampo += -2,"N",j-1, "N",j+1;      
        ampo += 1, s,j;
        ampo += 2*sign, s, j, "N",j+1;
        ampo += -2*sign, s, j, "N", j-1;
        auto Gi2 = toMPO(ampo);
        
        // Field operators
        ampo = AutoMPO(sites);
        ampo += 1,"N",j-1;
        ampo += -1*NE, "Id", j-1;
        auto Efield = toMPO(ampo);
        ampo = AutoMPO(sites);
        ampo += 1,"N",j-1,"N",j-1;
        ampo += NE*NE,"Id",j-1;
        ampo += -2*NE,"N",j-1;
        auto Efield2 = toMPO(ampo);
        
        //auto bondket = psi(j-1)*psi(j)*psi(j+1);
        Real gauge0 = inner(psiInit0,Gi,psiInit0);
        Real gauge = inner(psi0, Gi, psi0);
        Real gaugevar =  inner(psi0, Gi2, psi0)-gauge*gauge;
        Real E =  inner(psi0, Efield, psi0);
        Real E2 =  inner(psi0, Efield2, psi0)-E*E;
        printfln("%d %.12f %.12f %0.12f %0.12f %0.12f",j,gauge0, gauge, gaugevar, E, E2);
      }
  }
  return 0;
}

// int main(int argc, char* argv[]) {
//   if(argc < 2) 
//     { 
//     printfln("Usage: %s inputfile",argv[0]); 
//     return 0; 
//     }
//   schwinger(argv[1]);
//   return 0;
// }


