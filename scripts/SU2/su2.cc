#include "itensor/all.h"

using namespace itensor;

int su2(const char *inputfile, Args const& args = Args::global()){
  printfln("");
  auto input = InputGroup(inputfile,"input");
  auto N = input.getInt("N", 30);
  auto a = input.getReal("a", 0.1);
  auto m = input.getReal("m",1);
  auto g = input.getReal("g", 1);
  auto nstates = input.getInt("nstates", 1);
  auto weight  = input.getReal("weight", 100);

  Real t1 = 1;
  Real t2 = 1;
  Real t3 = 1;

  auto sites = Electron(N, {"ConserveQNs=",true, "ConserveSz",true});
  auto ampo = AutoMPO(sites);

  // 1d lattice with N spinful sites (i.e. N staggered sites, N/2 physical sites), labeled 1,...,N.

  // Interaction term H_I
  for(int j=1; j < N; j+=1)
    {
      ampo += t1*0.5/a, "Cdagup", j, "Cup", j+1;
      ampo += t1*0.5/a, "Cdagup", j+1, "Cup", j;
      ampo += t1*0.5/a, "Cdagdn", j, "Cdn", j+1;
      ampo += t1*0.5/a, "Cdagdn", j+1, "Cdn", j;
    }

  // Matter term H_M
  for(int j=1; j < N+1; j+=1)
    {
      ampo += t2*m*(2*(j%2)-1), "Cdagup", j, "Cup", j;
      ampo += t2*m*(2*(j%2)-1), "Cdagdn", j, "Cdn", j;
    }

  // Field term H_E
  for(int j=1; j < N+1; j+=1)
    {
    for(int k=j; k < N+1; k+=1) 
      {
	ampo += -t3*g*g*a/8*(k-j), "Cdagup", j, "Cup", j, "Cdagup", k, "Cup", k;
	ampo += -t3*g*g*a/8*(j-k), "Cdagup", j, "Cup", j, "Cdagdn", k, "Cdn", k;
	ampo += -t3*g*g*a/8*(j-k), "Cdagdn", j, "Cdn", j, "Cdagup", k, "Cup", k;
	ampo += -t3*g*g*a/8*(k-j), "Cdagdn", j, "Cdn", j, "Cdagdn", k, "Cdn", k;
	ampo += -t3*g*g*a/4*(k-j), "Cdagup", j, "Cdn", j, "Cdagdn", k, "Cup", k;
	ampo += -t3*g*g*a/4*(k-j), "Cdagdn", j, "Cup", j, "Cdagup", k, "Cdn", k;
      }
    }
  
  auto H = toMPO(ampo);

  auto sweeps     = Sweeps(15);
  sweeps.noise()  = 1E-6,1E-8,1E-10,1E-10,1E-10,1E-12,0,0,0,0,0,0;
  sweeps.niter() = 2;
  sweeps.maxdim() = 10,20,40,40,80,80,160,160,320,320,640,640,640,640,640;
  sweeps.cutoff() = 1E-10;
  Real  dE      = 1E-3;


  auto state = InitState(sites);
  for(int j=1; j < N+1; j++)
    {
      if(j%2)
	{
	  state.set(j,"Emp");
	}
      else
	{
	  state.set(j,"UpDn");
	}
    }
  auto psiInit0 = MPS(state);

  // Schmidt spectrum at half-cut
  // auto b = N-1;
  // psi0.position(b); 
  // auto r = leftLinkIndex(psi0,b);
  // auto t = siteIndex(psi0,b);
  // auto [U,S,V] = svd(psi0(b),{r,t});
  // auto u = commonIndex(U,S);
  // printfln("Shmidt coefficients:");
  // for(auto n : range1(dim(u)))
  //   {
  //     printfln("%.12f", elt(S,n,n));
  //   }

  // Excited states
  std::vector<MPS> states;
  std::vector<double> energies;

  for (int j = 0; j < nstates; j++) {
    auto [energy,psi] = dmrg(H, states, randomMPS(state), sweeps, {"Silent=", true, "Weight=", weight, "EnergyErrgoal=", dE});
    states.push_back(psi);
    energies.push_back(energy);
  }
  
  for (int j = 0; j < nstates; j++) {
    printfln("E%d= %d %.12f\n", j, energies[j]);
  }
  
  for(int k=0; k < 2; k+=1)
    {
      for(int j=1; j < N+1; j+=1)
	{
	  ampo = AutoMPO(sites);
	  ampo += 1, "Ntot", j;
	  auto nj   = toMPO(ampo);
	  Real occupation = inner(states[k], nj, states[k]);
	  printf("%.12f, ", occupation);  
	}
      printf("\n");
    }

  // Time evolution under the Hamiltonian
  // auto T = 0.1;
  // auto expH = toExpH(ampo,T*Cplx_i);
  // auto psi_T = applyMPO(expH,psiInit0,args);  
  return 0; 
}
