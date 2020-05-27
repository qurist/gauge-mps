#include "itensor/all.h"
#include "ladder.h"

using namespace itensor;

using Schwinger = MixedSiteSet<LadderSite, SpinHalfSite>;

int main(int argc, char* argv[]){
  if(argc < 2) 
    { 
    printfln("Usage: %s inputfile",argv[0]); 
    return 0; 
    }
  auto input = InputGroup(argv[1],"input");
  auto N = input.getInt("N", 30);
  auto NE = input.getInt("NE", 5);

  auto sites = Schwinger(2*N+1, {"ConserveSz=",true, "ConserveNb=", false, "MaxOcc=",2*NE});//SiteSet(N,4*NE+2);
  auto ampo = AutoMPO(sites);

  auto J = 1;
  auto mu = 5*J;
  auto x = 100*J;
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
      //ampo += -4*NE*l, "N", j-1;
      //ampo += -4*NE*l, "N", j+1;
      //ampo += 4*NE*NE*l, "Id", j+1;
      ampo += l, s,j;
      ampo += 2*sign*l, s, j, "N",j+1;
      ampo += -2*sign*l, s, j, "N", j-1;
    }
  
  auto H = toMPO(ampo);

  auto sweeps = Sweeps(10);
  sweeps.noise() = 1E-6,1E-8,1E-10,1E-12;
  sweeps.maxdim() = 10,20,50,100,200;
  sweeps.cutoff() = 1E-10;


  auto state = InitState(sites);
  for(auto j : range1(2*N+1))
    {
      if(j%4==2) state.set(j,"Dn");
      else if(j%4==0) state.set(j,"Up");
      else state.set(j,str(NE));
    }
  auto psiInit0 = MPS(state);

  auto [energy,psi0] = dmrg(H, psiInit0, sweeps, {"Quiet=",true});
  printfln("Ground State Energy = %.12f",energy/(2*N*x));

  //Excited state
  auto wfs = std::vector<MPS>(1);
  wfs.at(0) = psi0;

  // Here the Weight option sets the energy penalty for
  // psi1 having any overlap with psi0
  //
  // sweeps = Sweeps(5);
  // sweeps.noise() = 1E-6,1E-8,1E-10,1E-12;
  // sweeps.maxdim() = 10,20,40,50,50;
  // sweeps.cutoff() = 1E-10;

  // auto psiInit1 = randomMPS(state);
  // auto [en1,psi1] = dmrg(H,wfs,psiInit1,sweeps,{"Quiet=",true,"Weight=",20.0});
  // printfln("\nExcited State Energy = %.10f",en1);


  //  Check gauge conditions in the ground state
    // for(int j=3; j <= 2*N-3; j+=2)
    // {
    //   psi.position(j-1);
    //   auto Gi = -op(sites,"N",j-1)*op(sites,"Id",j)*op(sites,"Id",j+1) + op(sites,"Id",j-1)*op(sites,"Id",j)*op(sites,"N",j+1);
    //   if (j%4==1) {
    // 	Gi += - op(sites,"Id",j-1)*op(sites,"projUp",j)*op(sites,"Id",j+1);
    //   }
    //   else if (j%4==3) {
    // 	Gi += - op(sites,"Id",j-1)*op(sites,"projDn",j)*op(sites,"Id",j+1);
    //   }
    //   //auto Gi2 = Gi*Gi;
    //   auto bondket = psi(j-1)*psi(j)*psi(j+1);
    //   auto bondbra = dag(prime(bondket, "Ladder,Site,Ladder"));
    //   Real gauge = elt(bondbra*prime(Gi)*bondket);
    //   //Real gauge2 = elt(bondbra*Gi2*bondket);
    //   printfln("%d %.12f %0.12f",j,gauge,gauge);
    // }

  // // Construct Gauge operator G_1
  // ampo = AutoMPO(sites);
  // ampo += 1,"N",2;
  // ampo += -NE,"Id",2;
  // ampo += -1, "projUp",1; 
  // auto G1 = toMPO(ampo);
  
  // // Construct Gauge-squared G_1^2
  // ampo = AutoMPO(sites);
  // ampo += 1,"N",2, "N",2;
  // ampo += -2*NE, "N", 2;
  // ampo += NE*NE, "Id", 2;
  
  // ampo += 2*NE+1, "projUp",1;
  // ampo += -2, "projUp", 1, "N",2;
  // auto G12 = toMPO(ampo);
  // Real gauge1 = inner(psi, G1, psi);
  // Real gauge10 = inner(psi0,G1,psi0);
  // Real gauge12 =  inner(psi, G12, psi);
  // auto var1 = gauge12 - gauge1*gauge1;
  // printfln("%d %.12f %.12f %0.12f %0.12f %0.12f",1, gauge10, gauge1,var1);  
 
  
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
      //ampo += -4*NE*l, "N", j-1;
      //ampo += -4*NE*l, "N", j+1;
      //ampo += 4*NE*NE*l, "Id", j+1;

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

  
  // for (int i : range1(N)){
  //   std::vector<QN> QNs = {}; //[two's]
  //   auto qints = Index::qnstorage(1+maxOcc);
  //   for (int j : range(2*NE+1)) {
  //     qint[j] = QN({"E="+str(j-NE),j});
  //     QNs.push_back(qn);
  //   }
  //   sites(i) = Index(QNs, "Site,n="+str(i));
  // }
  // printFln(sites(3));
  return 0;
}
