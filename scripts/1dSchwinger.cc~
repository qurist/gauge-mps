#include "itensor/all.h"
#include "ladder.h"

using namespace itensor;

using Schwinger = MixedSiteSet<SpinHalfSite, LadderSite>;

int main(){
  int N = 10;
  int NE = 7;
  // N sites each equipped with a matter field (spin 1/2) and
  // electric field with cutoff NE (2*NE+1 dimensional ladder)
  // When symmetry is added, use a QN index set IS as follows:
  //   auto sites = SiteSet(IndexSet IS)

  auto sites = Schwinger(2*N-1, {"ConserveSz=",true, "ConserveNb=", false, "MaxOcc=",2*NE});//SiteSet(N,4*NE+2);

  auto ampo = AutoMPO(sites);
  
  auto J = 1;
  auto mu = 5;
  auto x = 100;
  auto l = 0;//10000;
  
  for(int j=1; j <= 2*N-3; j+=2)
    {
      // Pure Gauge Term
      ampo += J, "N",j+1,"N",j+1;
      ampo += -2*NE*J, "N",  j+1;
      ampo += NE*NE*J, "Id", j+1;
      // Matter Term
      auto s = "";
      if (j%4==1) {
	s = "projUp";
	  }
      else
	if (j%4==3) {
	  s = "projDn";
	}
      ampo += mu,s,j;
      // Matter-Gauge Term
      ampo += x,"S+",j,"N+",j+1,"S-",j+2;
      ampo += x,"S+",j+2,"N-",j+1,"S-",j;
      //ampo += x,"S+",j,"S-",j+2;
      //ampo += x,"S+",j+2,"S-",j;

      //Gauge Lagrangian sum_i lambda * [ Gi^2 ]
      if (j!=1){
      ampo += l,"N",j-1, "N",j-1;
      ampo += l,"N",j+1, "N",j+1;
      ampo += -2*l, "N", j-1, "N", j+1;
      ampo += l, s,j;
      ampo += -2*l, s, j, "N",j+1;
      ampo += 2*l, s, j, "N", j-1;
      }
    }
  
  auto H = toMPO(ampo);

  auto sweeps = Sweeps(10);
  //very important to use noise for this model
  sweeps.noise() = 1E-8,1E-8,1E-12,1E-15;
  sweeps.maxdim() = 10,20,50,100,200;
  sweeps.cutoff() = 1E-10;


  auto state = InitState(sites);
  for(auto j : range1(2*N-1))
    {
      if(j%4==1) state.set(j,"Dn");
      else if(j%4==3) state.set(j,"Up");
      else state.set(j,str(NE));
    }
  auto psi0 = MPS(state);

  auto [energy,psi] = dmrg(H, psi0, sweeps, {"Quiet=",true});
  printfln("Ground State Energy = %.12f",energy/(N*x));

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

  // Construct Gauge operator G_1
  ampo = AutoMPO(sites);
  ampo += 1,"N",2;
  ampo += -NE,"Id",2;
  ampo += -1, "projUp",1; 
  auto G1 = toMPO(ampo);
  
  // Construct Gauge-squared G_1^2
  ampo = AutoMPO(sites);
  ampo += 1,"N",2, "N",2;
  ampo += -2*NE, "N", 2;
  ampo += NE*NE, "Id", 2;
  
  ampo += 2*NE+1, "projUp",1;
  ampo += -2, "projUp", 1, "N",2;
  auto G12 = toMPO(ampo);
  Real gauge1 = inner(psi, G1, psi);
  Real gauge10 = inner(psi0,G1,psi0);
  Real gauge12 =  inner(psi, G12, psi);
  auto var1 = gauge12 - gauge1*gauge1;
  printfln("%d %.12f %.12f %0.12f %0.12f %0.12f",1, gauge10, gauge1,var1);  
 
  
  for(int j=3; j <= 2*N-3; j+=2)
    {
      psi.position(j-1);
      // Construct Gauge operator G_i
      ampo = AutoMPO(sites);
      ampo += -1,"N",j-1;
      ampo += 1,"N",j+1;
      if (j%4==1) {
    	ampo += -1, "projUp",j; 
      }
      else if (j%4==3) {
    	ampo+= -1, "projDn",j; 
      }
      auto Gi = toMPO(ampo);
      
      // Construct Gauge-squared G_i^2
      ampo = AutoMPO(sites);
      ampo += 1,"N",j-1, "N",j-1;
      ampo += 1,"N",j+1, "N",j+1;
      ampo += -2, "N", j-1, "N", j+1;
      auto s = "";
      if (j%4==1) {
	s = "projUp";
      }
      else if (j%4==3) {
	s = "projDn";
      }
      ampo += 1, s,j;
      ampo += -2, s, j, "N",j+1;
      ampo += 2, s, j, "N", j-1;
      
      auto Gi2 = toMPO(ampo);
      // Field operators
      ampo = AutoMPO(sites);
      ampo += 1,"N",j-1;
      auto Efield = toMPO(ampo);
      ampo = AutoMPO(sites);
      ampo += 1,"N",j-1,"N",j-1;
      auto Efield2 = toMPO(ampo);
      
      //auto bondket = psi(j-1)*psi(j)*psi(j+1);
      Real gauge = inner(psi, Gi, psi);
      Real gauge0 = inner(psi0,Gi,psi0);
      Real gauge2 =  inner(psi, Gi2, psi);
      Real E =  inner(psi, Efield, psi);
      Real E2 =  inner(psi, Efield2, psi)-E*E;
      auto var = gauge2 - gauge*gauge;
      printfln("%d %.12f %.12f %0.12f %0.12f %0.12f",j,gauge0, gauge,var, E, E2);
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
