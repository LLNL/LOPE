/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    scalarTransportFoam

Description
    Solves the steady or transient transport equation for a passive scalar.

\*---------------------------------------------------------------------------*/
#include <sys/stat.h>
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"
// #include "BridgeProblem.hpp"
// #include "SimplestNLP.hpp"
// #include "bridge.hpp"
#include "PstreamGlobals.H"
#include "IFstream.H"
//// #include "mfem/linalg/vector.hpp"
#include <map>
#include <vector>
#include <algorithm>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


inline bool file_exists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

struct dataTransfer
{   

    volScalarField * rFiberCell;     
    volScalarField * invPermeability;
    volScalarField * porosity;
    volScalarField * areaPerVol;
    volScalarField * dInvPermeability;
    volScalarField * dPorosity;
    volScalarField * dAreaPerVol;
    volScalarField * gamma;
    volScalarField * gammaAux;
    volScalarField * gammaRaw;        
    const cellList *  cells;
    const dimensionedScalar * rmin;
    const dimensionedScalar * rmax;
    dimensionedScalar * unitCellLen;
    List<double> * radiusInterpIn;
    List<double> * permInterpIn;    
    List<double> * porosityInterpIn;
    List<double> * areaInterpIn;
    List<double> * DpermInterpIn;    
    List<double> * DporosityInterpIn;
    List<double> * DareaInterpIn;
    const double * deltaRInterp;
    const volVectorField * centroids;      
    volScalarField * permeability;
    volScalarField * effDCO;
    volScalarField * effDCR;
    volScalarField * effKappa1;
    volScalarField * effKappa2;
    dimensionedScalar * permeabilityRef;
    dimensionedScalar * kappa1;    
    dimensionedScalar * kappa2;    
    dimensionedScalar * DCO;    
    dimensionedScalar * DCR;
    volVectorField * U;
    surfaceScalarField * phi;
    volScalarField * p;
    autoPtr<incompressible::turbulenceModel> * turbulence;
    label * pRefCell;
    scalar * pRefValue;
    dimensionedScalar * nu;
    fv::options * fvOptions;
    volVectorField * Ua;
    surfaceScalarField * phia;
    volScalarField * pa;
    singlePhaseTransportModel * laminarTransport;
    volTensorField * gradV;
    label * paRefCell;
    scalar * paRefValue;
    volScalarField * sens;                
    volScalarField * sensRaw;
    dimensionedScalar * rFilter;
    volScalarField * sensVel;
    double * betaTanh;
    dimensionedScalar * etaTanh;
    dimensionedScalar * min_gamma;                                 
    dimensionedScalar * max_gamma;
    const fvPatch * inletPatch;
    const fvPatch * topPatch;   
    double * brOptimizationObjective;
    double * brGammaStart;
    int * brGammaSize;    
    double * brGammaMin;
    double * brGammaMax;
    double * brGamma;
    double * brSensitivities;
    double * brLowerBound;
    double * brUpperBound;
    volScalarField * bR;
    volScalarField * bO;
    volScalarField * gR;
    volScalarField * gO;
    const dimensionedScalar * c1_km;
    const dimensionedScalar * c2_km;
    const dimensionedScalar * c3_km;
    const dimensionedScalar * c4_km;
    const dimensionedScalar * velRef;
    const dimensionedScalar * lenRef;
    volScalarField * CO;
    volScalarField * CR;
    volScalarField * Phi1;
    volScalarField * Phi2;
    volScalarField * COa;
    volScalarField * CRa;
    volScalarField * Psi1;
    volScalarField * Psi2;
    volScalarField * kmR;
    volScalarField * kmO;
    dimensionedScalar * alphaA;
    dimensionedScalar * alphaC;
    dimensionedScalar * CRef;
    dimensionedScalar * i0;
    dimensionedScalar * FRT;
    dimensionedScalar * nelec;
    dimensionedScalar * stdPotential;   
    dimensionedScalar * faradayConstant;
    volScalarField * deltaPhi;
    volVectorField * I1;
    volVectorField * I2;    
    dimensionedScalar * correctDim1UaEqn;
    dimensionedScalar * correctDim2UaEqn;
    volVectorField * velDerivInTerm;
    volVectorField * cOaGradCO;
    volVectorField * cRaGradCR;
    volScalarField * exchangeCurrent;
    volScalarField * sensConc;
    volScalarField * sensPhi;    
    volScalarField * dkmO_drf;
    volScalarField * dkmR_drf;
    volScalarField * din_dkmR;
    volScalarField * din_dkmO;
    volScalarField * din_drf;
    const volScalarField::Internal * meshVolumes;
    Foam::Time * runTime;
    double * VC;
    double * VCsens;
    double * VT;
}; 

struct dataTransfer transferVars;


void setNumberDecisionVars();

void setStartingValue();        

void setBounds(std::map<int,int> & map);        

//void saveResults(std::function<void()> &);

double evalObjective();

double evalConstraint();

void evalConstraintSensitivity();

void newDesignParameters();

void evalSensitivities();

/*

class SimpleProblem : public topopt::BridgeProblem
{
protected:
  int myrank_;
  int local_d_size_;
  MPI_Comm comm_;
  std::map<int,int> dv_to_field_; // maps indices in dv array to indices in gamma. Fixed indices are skipped over
public:
  // Adjust the local partition size during SetDecisionVariableBounds
  SimpleProblem(MPI_Comm comm, int myrank, std::map<int, int> &map) :
    myrank_(myrank),
    comm_(comm),
    local_d_size_(map.size()),
    dv_to_field_(map)
  {}
   
  virtual void SetDecisionVariableBounds(double * lower_bound, double * upper_bound)
  { 

    setBounds(dv_to_field_);
    std::transform(dv_to_field_.begin(),
		   dv_to_field_.end(),
		   lower_bound,
                   [&](const std::pair<const int,int> pair) { return (transferVars.brLowerBound)[pair.second]; });

    std::transform(dv_to_field_.begin(),
		   dv_to_field_.end(),
		   upper_bound,
                   [&](const std::pair<const int,int> pair) { return (transferVars.brUpperBound)[pair.second]; });

  }

   virtual void SetConstraintBounds(double * lower_bound, double * upper_bound)
   {
     lower_bound[0] = 0.0;
     upper_bound[0] = *transferVars.VT;
     std::cout << "JW: SetConstraints(" << myrank_ << "):" << std::endl;     
   }
      
  virtual void SetStartingPoint(double * starting_values)
  {
    setStartingValue();

    std::transform(dv_to_field_.begin(),
		   dv_to_field_.end(),
		   starting_values,
                   [&](const std::pair<const int,int> pair) { return (transferVars.brGammaStart)[pair.second]; });

    // std::copy(transferVars.brGammaStart,
    //           transferVars.brGammaStart+local_d_size_,
    //           starting_values);
    // std::cout << "\n ..: Print each des var \n";
    // for (int i = 0; i < 10; i++)
    // {
    //   std::cout << starting_values[i] << "\n";
    // }
  }
            
  virtual void EvaluateObjectiveGradients(double *dp) 
  {
    // evalSensitivities();  // this is done on all procs
    std::transform(dv_to_field_.begin(),
		   dv_to_field_.end(),
		   dp,
                   [&](const std::pair<const int,int> pair) { return (transferVars.brSensitivities)[pair.second]; });

    
    // std::copy(transferVars.brSensitivities,
    //           transferVars.brSensitivities+local_d_size_,
    //           dp);
    // std::cout << ".. dp gradients \n";
    // for (int i = 0; i < 100; i++)
    // {
    //   std::cout << "..: Grad: " << dp[i] << " " << transferVars.brSensitivities[i] << "\n";
    // }
    // std::cout << "JW: EvalObjectiveGradients(" << myrank_ << "):" 
    //           <<  transferVars.brSensitivities[0] << std::endl;
    
  }
      
   virtual double EvaluateObjective() 
   {
      double obj = evalObjective();  // this is done on all procs
      std::cout << "JW: EvaluateObjective(" << myrank_ << "):" << obj << std::endl;
      return obj;
   }

   virtual void EvaluateConstraints(double * g)
   {
     
     std::cout << "JW: EvalConstraints(" << myrank_ << "):" << std::endl;
     double my_g = evalConstraint();
     g[0] = my_g;


   }
   
   virtual void EvaluateConstraintGradients(std::vector<double *> &Gradient)
   {
     
    //  std::cout << "JW: EvalConstraintGradients(" << myrank_ << "):" << std::endl;
    evalConstraintSensitivity();

    std::transform(dv_to_field_.begin(),
		   dv_to_field_.end(),
		   Gradient[0],
                   [&](const std::pair<const int,int> pair) { return (transferVars.VCsens)[pair.second]; });

    
    // std::copy(transferVars.VCsens,
    //           transferVars.VCsens+local_d_size_,
    //           Gradient[0]);
   }

  mutable std::function<void()> dump_callback;
  
  virtual void NewXCallback()
  {
    std::vector<double> dv(local_d_size_);
    GetDecisionVariables(dv.data());
    //    transferVars.brGamma
    for (const auto &pair : dv_to_field_)
      {
	transferVars.brGamma[pair.second] = dv[pair.first];
      }
    
    newDesignParameters();  // this is done on all procs
    evalSensitivities();

  }
  
  virtual void IterateCallback()
  {
    
    Info << "Iterate Call Back." << endl;
    saveResults(dump_callback); 
  }


};
*/



int main(int argc, char *argv[])
{
           
                //#include "setRootCase.H"   
    Foam::argList args(argc, argv); //this constructor will also initialze MPI?
                                    //Constructor calls "runPar"... 
                                    //I think this happens at line 71 of parRun.H
    if (!args.checkRootCase())
    {
        Foam::FatalError.exit();
    } 

    
    //#include "createTime.H"
    Foam::Info<< "Create time\n" << Foam::endl;
    Foam::Time runTime(Foam::Time::controlDictName, args);

    // Info << runTime << endl;

    
    
    //#include "createMesh.H"
    Foam::Info
    << "Create mesh for time = "
    << runTime.timeName() << Foam::nl << Foam::endl;

    Foam::fvMesh mesh
    (
        Foam::IOobject
        (
            Foam::fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );
    
    simpleControl simple(mesh);

    #include "createFields.H"
    

    fv::options& fvOptions(fv::options::New(mesh));

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    
    dimensionedScalar correctDim1UaEqn("correctDim1UaEqn", 
                dimensionSet(0, -1, -3, 0, 0, 0, 0), 1.0);

    dimensionedScalar correctDim2UaEqn("correctDim2UaEqn", 
                dimensionSet(0, -2, -2, 0, 0, 0, 0), 1.0);
    
    dimensionedScalar correctDim1("correctDim1", 
                dimensionSet(-1, 1, 3, 0, 0, 0, 0),1.0);

    dimensionedScalar correctDim2("correctDim2", 
                dimensionSet(0, -2, 3, 0, 0, 0, 0),1.0);

    dimensionedScalar correctDim3("correctDim3", 
                dimensionSet(-1, -3, 3, 0, 0, 2, 0),1.0);

    dimensionedScalar correctDimVelMax("correctDimVelMax", 
                dimensionSet(0, -6, 5, 0, 0, 0, 0),1.0);
                
    dimensionedScalar min_gamma("min_gamma",
                                 dimless, 1.0e-4);
    dimensionedScalar max_gamma("max_gamma",
                                 dimless, 1.0 - 1.0e-4);

    dimensionedScalar min_gammaRaw("min_gammaRaw",
                                 dimless, 0.0);
    dimensionedScalar max_gammaRaw("max_gammaRaw",
                                 dimless, 1.0);

    volVectorField cOaGradCO(COa*correctDim1UaEqn*fvc::grad(CO));
    volVectorField cRaGradCR(CRa*correctDim1UaEqn*fvc::grad(CR));    
    volTensorField gradV(fvc::grad(U));     
    
    volVectorField velDerivInTerm
                   (                  
                      areaPerVol * nelec * faradayConstant
                    * (
                          CO * Foam::exp(-1.0*alphaC * FRT * deltaPhi) 
                        - CR * Foam::exp(alphaA * FRT * deltaPhi)
                      ) 
                    /
                      Foam::pow(
                          nelec * faradayConstant * CRef / i0
                        + Foam::exp(-1.0*alphaC * FRT * deltaPhi)
                          / (gO * Foam::pow(mag(U/velRef)+1.e-20,bO))
                        + Foam::exp(alphaA * FRT * deltaPhi)
                          / (gR * Foam::pow(mag(U/velRef)+1.e-20,bR))
                         ,2)
                    * (
                          Foam::exp(alphaA * FRT * deltaPhi)*bR
                          / (gR * Foam::pow(mag(U/velRef)+1.e-20,bR+1.))
                        + Foam::exp(-1.0*alphaC * FRT * deltaPhi)
                          / (gO * Foam::pow(mag(U/velRef)+1.e-20,bO+1.))
                      )
                    * ( (CRa-COa) / nelec / faradayConstant + Psi1 - Psi2 )
                    * correctDim2UaEqn  
                    * U/(mag(U) + 1.e-100*velRef)                   
                   );
    
      
    const label& patchID = mesh.boundaryMesh().findPatchID("top"); 
    const fvPatch& topPatch = mesh.boundary()[patchID];     
    const label& patchID2 = mesh.boundaryMesh().findPatchID("btm"); 
    const fvPatch& btmPatch = mesh.boundary()[patchID2]; 
    const label& patchID3 = mesh.boundaryMesh().findPatchID("outlet"); 
    const fvPatch& outletPatch = mesh.boundary()[patchID3];        
    const label& patchID4 = mesh.boundaryMesh().findPatchID("inlet"); 
    const fvPatch& inletPatch = mesh.boundary()[patchID4];        
               
    const scalar meshVol = gSum(mesh.V());
        
    double timeInterp = runTime.value();
        
    Info << "Initializing Interpolation Parameters at Time = " 
         << runTime.timeName() << endl;
    // brGammaSize
    // Initialize interpolation parameters
    
    
    double betaTanh = betaTanh0_.value();    
    double qArea = qArea0_.value();    
    double qPerm = qPermInterp0_.value();
    double gammaTrnPerm = gammaTrnPermInterp0_.value();   

    Info << "qArea: " << qArea << endl;        
    Info << "qPerm: " << qPerm << endl;
    Info << "betaTanh: " << betaTanh << endl;  
    Info << "gammaTrnPerm: " << gammaTrnPerm << endl;      
    
    // Variale step-length parameter
    
    double maxStepLength = stepLength.value();
    int restepFlag = 0;
    
    //stepLength = stepLengthInit;


    List<double> radiusInterpIn(IFstream("constant/radiusInterpIn.txt")());
    List<double> permInterpIn(IFstream("constant/permInterpIn.txt")());    
    List<double> porosityInterpIn(IFstream("constant/porosityInterpIn.txt")());
    List<double> areaInterpIn(IFstream("constant/areaInterpIn.txt")());
    List<double> DpermInterpIn(IFstream("constant/dPermInterpIn.txt")());    
    List<double> DporosityInterpIn(IFstream("constant/dPorosityInterpIn.txt")());
    List<double> DareaInterpIn(IFstream("constant/dAreaInterpIn.txt")());    
    

    const dimensionedScalar rmin = printResLen/2.0;
    const dimensionedScalar rmax = (unitCellLen - 2.0 * printResLen) / 4.0;
    const double deltaRInterp = radiusInterpIn[1] - radiusInterpIn[0];

    
    int valueN = 0;
    forAll(radiusInterpIn, valueI)
    {
        valueN = valueI;
    }
    
    
    Info << "rmin: " << rmin.value() << endl;
    Info << "rmax: " << rmax.value() << endl;
    Info << "radiusInterpIn entry 0: " << radiusInterpIn[0] << endl;
    Info << "radiusInterpIn entry " << valueN << ": " 
                                    << radiusInterpIn[valueN] << endl;

    if ( rmin  < radiusInterpIn[0]*unitCellLen
               || rmax > radiusInterpIn[valueN] * unitCellLen )
    {
        Info << "INTERPOLATINN ERROR! " 
             << "SIMULATION BOUNDS EXCEED INTERPOLATION DATA" << endl;
    }  
    
    
    volScalarField dkmO_drf = kmO/rFiber;
    volScalarField dkmR_drf = kmR/rFiber;
    volScalarField din_dkmR = i0/kmR;
    volScalarField din_dkmO = i0/kmO;
    volScalarField din_drf = exchangeCurrent/rFiber;


    ///Initialize transferVars
    double brOptimizationObjective = 1.0;
    int brGammaSize = gamma.size(); 
    double VC = 0.0;
    double VT = 1.1;
    transferVars =
    {   
        &rFiberCell,     
        &invPermeability,
        &porosity,
        &areaPerVol,    
        &dInvPermeability,
        &dPorosity,
        &dAreaPerVol,
        &gamma,
        &gammaAux,
        &gammaRaw,        
        &mesh.cells(),
        &rmin,
        &rmax,
        &unitCellLen,
        &radiusInterpIn,
        &permInterpIn,    
        &porosityInterpIn,
        &areaInterpIn,
        &DpermInterpIn,    
        &DporosityInterpIn,
        &DareaInterpIn,
        &deltaRInterp,
        &mesh.C(),
        &permeability,
        &effDCO,
        &effDCR,
        &effKappa1,
        &effKappa2,
        &permeabilityRef,
        &kappa1,    
        &kappa2,    
        &DCO,    
        &DCR,
        &U,
        &phi,
        &p,
        &turbulence,
        &pRefCell,
        &pRefValue,
        &nu,
        &fvOptions,
        &Ua,
        &phia,
        &pa,
        &laminarTransport,
        &gradV,
        &paRefCell,
        &paRefValue,
        &sens,                
        &sensRaw,
        &rFilter,
        &sensVel,
        &betaTanh,
        &etaTanh,
        &min_gamma,                                 
        &max_gamma,
        &inletPatch,
        &topPatch,
        &brOptimizationObjective,
        new double [brGammaSize],
        &brGammaSize,
        &min_gamma.value(), //min_gammaRaw.value(),
        &max_gamma.value(), //max_gammaRaw.value(),
        new double [brGammaSize],
        new double [brGammaSize],
        new double [brGammaSize],
        new double [brGammaSize],
        &bR,
        &bO,
        &gR,
        &gO,
        &c1_km,
        &c2_km,
        &c3_km,
        &c4_km,
        &velRef,
        &lenRef,
        &CO,
        &CR,
        &Phi1,
        &Phi2,
        &COa,
        &CRa,
        &Psi1,
        &Psi2,
        &kmR,
        &kmO,
        &alphaA,
        &alphaC,
        &CRef,
        &i0,
        &FRT,
        &nelec,
        &stdPotential,
        &faradayConstant,
        &deltaPhi,
        &I1,
        &I2,
        &correctDim1UaEqn,
        &correctDim2UaEqn,
        &velDerivInTerm,
        &cOaGradCO,
        &cRaGradCR,
        &exchangeCurrent,
        &sensConc,
        &sensPhi,
        &dkmO_drf,
        &dkmR_drf,
        &din_dkmR,
        &din_dkmO,
        &din_drf,
        &mesh.V(),
        &runTime,
        &VC,
        new double [brGammaSize],  // VCsens
        &VT
    };


    U.storePrevIter();
    p.storePrevIter();
   
    setNumberDecisionVars();
    std::cout << "  : Set Number Decision vars \n";
    std::map<int, int> map;
    setBounds(map);
    std::cout << "  : Set Bounds \n";
    setStartingValue();
    std::cout << "  : Set starting value \n";
    newDesignParameters();

    std::fill_n (transferVars.brSensitivities, brGammaSize, 0.0);

    evalSensitivities();
    double initial_value =  evalObjective();


    evalConstraint();
    evalConstraintSensitivity();
    
    //
    /*
    std::fill_n (transferVars.brSensitivities, brGammaSize, 0.0);
    
    Info << "Size of MPICommuicators: " << Foam::PstreamGlobals::MPICommunicators_.size() << endl;  
    Info << "MPICommiunicators: " << Foam::PstreamGlobals::MPICommunicators_ << endl;
    Info << "Pstream::masterNo: " << Pstream::masterNo() << endl;    
    
    
    MPI_Comm mpicomm;
    if (Foam::PstreamGlobals::MPICommunicators_.size() == 0) 
    {
        mpicomm = MPI_COMM_WORLD;
        MPI_Init(NULL, NULL);
    } else {
        mpicomm = Foam::PstreamGlobals::MPICommunicators_[0];
    }

    int my_rank = 0;
    MPI_Comm_rank( mpicomm, &my_rank);

    int local_ndv = brGammaSize;

    SimpleProblem *s =  new SimpleProblem(mpicomm, my_rank, map);
 
    std::cout << "  : SimpleProblem initialized! \n";
    double x[map.size()];
    std::transform(map.begin(),
      map.end(),
      x,
      [&](const std::pair<const int,int> pair) {
          return transferVars.brGamma[pair.second];
        });
    s->SetDecisionVariables(x, map.size());
    std::cout << "  : Decision Variables initialized! \n";
    s->SetNumberOfNLConstraints(1); 
    bridge_optimizer_settings settings;
    settings.acceptable_tol = 1e-30;
    settings.tol = 1e-30;
    settings.fabstol = 0;
    settings.max_it = 5000;
    settings.freltol = 0;
    settings.test_deriv = false;
    settings.dh = 1e-4;  // derivative step size 1e-2 gave good rel error?
    //settings.optimizer_solver = "ipopt";
    settings.optimizer_solver = "mma";
    // mma relative and abs tolerance options
    settings.numeric_options["f_abs_tol"] = 0.0;
    settings.numeric_options["f_rel_tol"] = 0.0;

    // Arbitrary ipopt options settings:
    // settings.integer_options["max_iter"] = 2;
    // settings.string_options["mu_strategy"] = "monotone";
    settings.string_options["derivative_test_print_all"] = "yes";
    // settings.string_options["accept_every_trial_step"] = "yes";
    // settings.string_options["nlp_scaling_method"] = "none";
    settings.string_options["nlp_scaling_method"] = "equilibration-based";
    // settings.numeric_options["obj_scaling_factor"] = 1e2;

    settings.string_options["linear_solver"] = "ma57";
    // settings.string_options["accept_every_trial_step"] = "yes";
    // settings.integer_options["accept_after_max_steps"] = 10;
    // settings.integer_options["limited_memory_max_history"] = 20;
    // settings.numeric_options["tiny_step_tol"] = 1e-5;
    // settings.numeric_options["tiny_step_y_tol"] = 1e-10;
    // settings.integer_options["watchdog_shortened_iter_trigger"] = 1;

    if (file_exists("MMA_Restart/state.Info.txt"))
      {
	// state file exists so restart
	settings.restart = true;
	Info << "Found restart file " << endl;
      }

    
    std::cout << "  : Call SimplestNLP wrapper \n";
    SimplestNLP wrapper(*s, mpicomm, 0);

    std::cout << "  : Start the optimization \n";
    SolveProblem(&wrapper, settings, s->dump_callback);


    s->GetDecisionVariables(x);
    for (int i = 0; i < 100 ; i++)
      Info << "Gamma Bulk Lido: " << x[i] << endl;
  
  
    delete [] transferVars.brGammaStart;
    delete [] transferVars.brGamma;
    delete [] transferVars.brSensitivities;
    delete [] transferVars.brUpperBound;
    delete [] transferVars.brLowerBound;
    delete [] transferVars.VCsens;
    delete s;
    return 0;
  */
  //
  
     while (simple.loop())
     {
         Info << nl << "Time = " << runTime.timeName() << nl << endl;                        
         timeInterp = runTime.value();
        

         forAll(mesh.C(), cellI)
         {
            transferVars.brGamma[cellI] -= transferVars.brSensitivities[cellI] * stepLength.value();
            if (mesh.C()[cellI][1] < 0.0)
            {
              transferVars.brGamma[cellI] = 1.0;
            
            }        
         }
    
                
         newDesignParameters();
         evalSensitivities();
         //brOptimizationObjective = evalObjective();                                              
        
             
         runTime.write();  
        
         const fvPatchField<vector>& topI2 =
                 topPatch.lookupPatchField<volVectorField, vector>("I2");                        
         Info << "Total Current at top: "
             << gSum( (topI2 & topI2.patch().Sf()) ) << endl;                   
        
         const fvPatchField<vector>& btmI1 =
                 btmPatch.lookupPatchField<volVectorField, vector>("I1");        
         Info << "Total Current at btm: "
             << gSum( (btmI1 & btmI1.patch().Sf()) ) << endl;    
         const fvPatchField<scalar>& topPhi2 = 
                 topPatch.lookupPatchField<volScalarField, scalar>("Phi2");
         const fvPatchField<scalar>& inletPressure =                 
                 inletPatch.lookupPatchField<volScalarField, scalar>("p");
         const fvPatchField<scalar>& outletCO = 
                     outletPatch.lookupPatchField<volScalarField, scalar>("CO");
         const fvPatchField<scalar>& outletCR = 
                     outletPatch.lookupPatchField<volScalarField, scalar>("CR");                
         const fvPatchField<scalar>& inletCO =   
                    inletPatch.lookupPatchField<volScalarField, scalar>("CO");
         const fvPatchField<scalar>& inletCR = 
                    inletPatch.lookupPatchField<volScalarField, scalar>("CR"); 
         const fvsPatchField<scalar>& inPhi = 
                    inletPatch.lookupPatchField<surfaceScalarField, scalar>("phi");
        
         Info << "Inlet flow (mol/s) of CO: " 
              << gSum(inletCO * phi.boundaryField()[patchID4]) << endl;
         Info << "Inlet flow (mol/s) of CR: " 
              << gSum(inletCR * phi.boundaryField()[patchID4]) << endl;                
         Info << "Outlet flow (mol/s) of CO: " 
              << gSum(outletCO * phi.boundaryField()[patchID3]) << endl;
         Info << "Outlet flow (mol/s) of CR: " 
              << gSum(outletCR * phi.boundaryField()[patchID3]) << endl;             
         Info << "Flow rate out (m^3/s): " 
              << gSum(phi.boundaryField()[patchID3]) << endl;
         Info << "Flow rate in (m^3/s): " 
              << gSum(phi.boundaryField()[patchID4]) << endl;                     
         Info << "Flow Power (J/s): " 
              << gSum(-1000.
                      *inletPressure*inPhi) << endl;
         Info << "Electric Power (J/s): " 
              << gSum((topI2 & topI2.patch().Sf())
                      *(-1.*topPhi2 - stdPotential.value()) ) << endl;        
         Info << "Total Power Loss (J/s): " 
              << gSum((topI2 & topI2.patch().Sf())
                      *(-1.*topPhi2 - stdPotential.value()) )
                 +gSum(-1000.
                      *inletPressure*inPhi) << endl;                 
    }
    Info<< "End\n" << endl;

    delete [] transferVars.brGammaStart;
    delete [] transferVars.brGamma;
    delete [] transferVars.brSensitivities;
    delete [] transferVars.brUpperBound;
    delete [] transferVars.brLowerBound;
    delete [] transferVars.VCsens;
    //delete s;
    return 0;

}
/*
void saveResults(std::function<void()> & dump)
{
    Foam::Time              & myRunTime         = *transferVars.runTime;
    volScalarField          & gamma             = *transferVars.gamma;
    volScalarField          & gammaRaw             = *transferVars.gammaRaw;
    double                  & brOptimizationObjective = *transferVars.brOptimizationObjective;       
    double                  & VC                = *transferVars.VC;       
    // Info << nl << "saveResults: Time = " << myRunTime.timeName() << nl << endl;
    ++myRunTime;
    Info << nl << "saveResults: Time New = " << myRunTime.timeName() << " objective: " << brOptimizationObjective << " constraint: " << VC <<endl;
    //gamma.write();  // only write gamma
    dump();
    // gammaRaw.write();
    myRunTime.write();
}
*/

double evalConstraint()
{

  volScalarField          & gamma             = *transferVars.gamma;
  const volScalarField::Internal & meshVolumes       = *transferVars.meshVolumes;
  double                  & VC                = *transferVars.VC;       

  VC = gSum(gamma*1.0*meshVolumes) / gSum(meshVolumes);
  std::cout << "evalConstraint: Constraint: " << VC << " \n";

  return VC;
}

void evalConstraintSensitivity()
{
  const volScalarField::Internal & meshVolumes       = *transferVars.meshVolumes;
  double                  & VC                = *transferVars.VC;
  double                  & VT                = *transferVars.VT;
  const cellList          & cells             = *transferVars.cells;

  double vc_vt = 1./gSum(meshVolumes);

  forAll(cells, cellI)
  {
    transferVars.VCsens[cellI] = vc_vt*meshVolumes[cellI]; // was 0.01?
  }

}

double evalObjective()
{

    surfaceScalarField      & phi               = *transferVars.phi;
    volScalarField          & p                 = *transferVars.p;
    const fvPatch           & inletPatch        = *transferVars.inletPatch;  
    const fvPatch           & topPatch          = *transferVars.topPatch;     
    double                  & brOptimizationObjective = *transferVars.brOptimizationObjective;       
    volScalarField          & Phi2              = *transferVars.Phi2;
    dimensionedScalar       & stdPotential      = *transferVars.stdPotential;
    Foam::Time              & myRunTime         = *transferVars.runTime;
    volScalarField          & gamma             = *transferVars.gamma;



    const fvPatchField<scalar>& inletPressure =                 
                inletPatch.lookupPatchField<volScalarField, scalar>("p");
    const fvsPatchField<scalar>& inPhi = 
                   inletPatch.lookupPatchField<surfaceScalarField, scalar>("phi");                   
    const fvPatchField<scalar>& topPhi2 = 
                topPatch.lookupPatchField<volScalarField, scalar>("Phi2");
    const fvPatchField<vector>& topI2 =
                topPatch.lookupPatchField<volVectorField, vector>("I2");  

        
    //brOptimizationObjective = gSum(-1000.*inletPressure*inPhi);
    
    brOptimizationObjective = gSum((topI2 & topI2.patch().Sf())
                        *(-1.*topPhi2 - stdPotential.value()) )
                        + gSum(-1000.*inletPressure*inPhi);
    //Info << "Flow Power (J/s): " << brOptimizationObjective << endl;
    std::cout << "BrObjective: " << brOptimizationObjective << " \n";
    return brOptimizationObjective;
}

void setNumberDecisionVars()
{
    int                     & brGammaSize       = *transferVars.brGammaSize;    
    volScalarField          & gamma             = *transferVars.gamma;
    
    
    brGammaSize = gamma.size();

}


void setStartingValue()
{
    volScalarField          & gammaRaw          = *transferVars.gammaRaw;
    const cellList          & cells             = *transferVars.cells;
    const volVectorField    & centroids         = *transferVars.centroids; 

    /// Initialize brGammaStart and brGamma
    //forAll(cells, cellI)
    //{
    // for (int i = 0; i < *transferVars.brGammaSize; i++){
	  //     transferVars.brGammaStart[i] = 0.6;
    //     transferVars.brGamma[i] = 0.6;        
    // }
    forAll(cells, cellI)
    {
        transferVars.brGammaStart[cellI] = 0.6;
        transferVars.brGamma[cellI] = 0.6;
        if (centroids[cellI][1] < 0.0)
        {
            transferVars.brGammaStart[cellI] = 1.0;
            transferVars.brGamma[cellI] = 1.0;
            // set lower bound equal to upper bound in this case
        }        
    }
    
       
}

void setBounds(std::map<int,int> & map)
{
    const cellList          & cells             = *transferVars.cells;
    dimensionedScalar       & min_gamma         = *transferVars.min_gamma;                                 
    dimensionedScalar       & max_gamma         = *transferVars.max_gamma;
    double                  & brGammaMin        = *transferVars.brGammaMin;
    double                  & brGammaMax        = *transferVars.brGammaMax;
    // double                  & lowerBound        = *transferVars.brLowerBound;
    // double                  & upperBound        = *transferVars.brUpperBound;
    const volVectorField    & centroids         = *transferVars.centroids; 

    std::vector<int> map_mod(*transferVars.brGammaSize, 0);

    //brGammaMin = min_gamma.value();
    //brGammaMax = max_gamma.value();
    
    brGammaMin = 0.0;
    brGammaMax = 1.0;
    /// First, update design variables with output of optimizer
    forAll(cells, cellI)
    {
        transferVars.brLowerBound[cellI] = brGammaMin;
        transferVars.brUpperBound[cellI] = brGammaMax;
        if (centroids[cellI][1] < 0.0)
        {
            transferVars.brLowerBound[cellI] = brGammaMax;
            map_mod[cellI] = 1;
            // set lower bound equal to upper bound in this case
        }        
    }
    // adjust local indicies so only non-fixed decision variables are sent to the optimizer
    map.clear();
    int counter = 0;
    for (int i = 0; i < map_mod.size(); i++)
      {
        if (map_mod[i] == 0)
          map[counter++] = i;
      }
}


void evalSensitivities()
{
    volVectorField          & U                 = *transferVars.U; 
    volVectorField          & Ua                = *transferVars.Ua;
    dimensionedScalar       & nu                = *transferVars.nu;
    volScalarField          & dInvPermeability  = *transferVars.dInvPermeability;    
    const dimensionedScalar & rmin              = *transferVars.rmin;
    const dimensionedScalar & rmax              = *transferVars.rmax;
    volScalarField          & gamma             = *transferVars.gamma;
    volScalarField          & sens              = *transferVars.sens;                
    volScalarField          & sensRaw           = *transferVars.sensRaw;
    dimensionedScalar       & rFilter           = *transferVars.rFilter;
    volScalarField          & sensVel           = *transferVars.sensVel;
    double                  & betaTanh          = *transferVars.betaTanh;
    dimensionedScalar       & etaTanh           = *transferVars.etaTanh;
    dimensionedScalar       & min_gamma         = *transferVars.min_gamma;                                 
    dimensionedScalar       & max_gamma         = *transferVars.max_gamma;
    const cellList          & cells             = *transferVars.cells;
    volScalarField          & exchangeCurrent   = *transferVars.exchangeCurrent;
    dimensionedScalar       & i0                = *transferVars.i0;
    volScalarField          & CO                = *transferVars.CO;
    volScalarField          & CR                = *transferVars.CR;
    volScalarField          & Phi1              = *transferVars.Phi1;
    volScalarField          & Phi2              = *transferVars.Phi2;
    volScalarField          & COa               = *transferVars.COa;
    volScalarField          & CRa               = *transferVars.CRa;
    volScalarField          & Psi1              = *transferVars.Psi1;
    volScalarField          & Psi2              = *transferVars.Psi2;
    dimensionedScalar       & alphaA            = *transferVars.alphaA;
    dimensionedScalar       & alphaC            = *transferVars.alphaC;
    dimensionedScalar       & FRT               = *transferVars.FRT;
    volScalarField          & deltaPhi          = *transferVars.deltaPhi;
    dimensionedScalar       & nelec             = *transferVars.nelec;
    dimensionedScalar       & faradayConstant   = *transferVars.faradayConstant;
    volScalarField          & kmR               = *transferVars.kmR;
    volScalarField          & kmO               = *transferVars.kmO;    
    dimensionedScalar       & CRef              = *transferVars.CRef;
    volScalarField          & rFiberCell        = *transferVars.rFiberCell;
    volScalarField          & dPorosity         = *transferVars.dPorosity;
    volScalarField          & dAreaPerVol       = *transferVars.dAreaPerVol;
    dimensionedScalar       & kappa1            = *transferVars.kappa1;
    dimensionedScalar       & kappa2            = *transferVars.kappa2;
    dimensionedScalar       & DCO               = *transferVars.DCO;
    dimensionedScalar       & DCR               = *transferVars.DCR;
    volScalarField          & porosity          = *transferVars.porosity;
    volScalarField          & areaPerVol        = *transferVars.areaPerVol;   
    volScalarField          & sensConc          = *transferVars.sensConc;  
    volScalarField          & sensPhi           = *transferVars.sensPhi;  
    const dimensionedScalar & lenRef            = *transferVars.lenRef;
    const dimensionedScalar & c2_km             = *transferVars.c2_km;
    const dimensionedScalar & c3_km             = *transferVars.c3_km;
    const dimensionedScalar & c4_km             = *transferVars.c4_km;
    volScalarField          & dkmO_drf          = *transferVars.dkmO_drf;
    volScalarField          & dkmR_drf          = *transferVars.dkmR_drf;
    volScalarField          & din_dkmR          = *transferVars.din_dkmR;
    volScalarField          & din_dkmO          = *transferVars.din_dkmO;
    volScalarField          & din_drf           = *transferVars.din_drf;
    const volVectorField    & centroids         = *transferVars.centroids; 
    scalar                  sumSens; 
    const volScalarField::Internal          & meshVol           = *transferVars.meshVolumes;
    volScalarField          & gammaAux          = *transferVars.gammaAux;
    
     
    dimensionedScalar correctDim2("correctDim2",
                                  dimensionSet(0, -2, 3, 0, 0, 0, 0),1.0);        
        
        exchangeCurrent =  i0 *  
                  (
                     CR * Foam::exp(alphaA * FRT * deltaPhi)
                   - CO * Foam::exp(-1.0*alphaC * FRT * deltaPhi)
                  )
                  /
                  ( 1. + i0/nelec/faradayConstant/CRef
                       * (
                             Foam::exp(alphaA * FRT * deltaPhi ) / kmR
                           + Foam::exp(-1.0*alphaC * FRT * deltaPhi) / kmO
                         )
                  );


        dkmO_drf =  kmO* (  -1.0/(rFiberCell+1.0e-20*lenRef) 
                            + (
                                  c4_km * Foam::log(mag(U)*rFiberCell/DCO+1.0e-20) 
                                - c2_km/(porosity+1.0e-20)
                              )*dPorosity
                         + (c3_km +c4_km*porosity)/(rFiberCell+1.0e-20*lenRef)
                         );
      

        dkmR_drf =  kmR* (  -1.0/(rFiberCell+1.0e-20*lenRef) 
                            + (
                                  c4_km * Foam::log(mag(U)*rFiberCell/DCR+1.0e-20) 
                                - c2_km/(porosity+1.0e-20)
                              )*dPorosity
                         + (c3_km +c4_km*porosity)/(rFiberCell+1.0e-20*lenRef)
                         );

        din_dkmO = CRef * Foam::pow(i0/CRef/kmO,2) 
                        * Foam::exp(-1.0*alphaC * FRT * deltaPhi)
                          /nelec/faradayConstant
                        * 
                          (
                              CR * Foam::exp(alphaA * FRT * deltaPhi)
                            - CO * Foam::exp(-1.0*alphaC * FRT * deltaPhi)
                          )
                        /
                        Foam::pow( 1. + i0/nelec/faradayConstant/CRef
                       * (
                             Foam::exp(alphaA * FRT * deltaPhi ) / kmR
                           + Foam::exp(-1.0*alphaC * FRT * deltaPhi) / kmO
                         ),2);

        din_dkmR = CRef * Foam::pow(i0/CRef/kmO,2) 
                        * Foam::exp(alphaA * FRT * deltaPhi)
                          /nelec/faradayConstant
                        * 
                          (
                              CR * Foam::exp(alphaA * FRT * deltaPhi)
                            - CO * Foam::exp(-1.0*alphaC * FRT * deltaPhi)
                          )
                        /
                        Foam::pow( 1. + i0/nelec/faradayConstant/CRef
                       * (
                             Foam::exp(alphaA * FRT * deltaPhi ) / kmR
                           + Foam::exp(-1.0*alphaC * FRT * deltaPhi) / kmO
                         ),2);

        din_drf = din_dkmO * dkmO_drf + din_dkmR * dkmR_drf;

        sensPhi =   (rmin - rmax) 
                    * (
                          dAreaPerVol * exchangeCurrent 
                        + areaPerVol * din_drf
                      )                  
                    * (Psi1 - Psi2)                    
                    - 1.5 * (rmin - rmax) 
                          * sqrt( 1.0 - porosity ) * dPorosity * kappa1
                          * (fvc::grad(Phi1)&fvc::grad(Psi1))
                    + 1.5 * (rmin - rmax) 
                          * sqrt( porosity ) * dPorosity * kappa2                    
                          * (fvc::grad(Phi2)&fvc::grad(Psi2));
                                    
        sensConc =  (rmin - rmax)
                    * (
                          dAreaPerVol * exchangeCurrent 
                        + areaPerVol * din_drf
                      )
                    * (CRa - COa)/nelec/faradayConstant
                    + 1.5 * (rmin - rmax) 
                          * sqrt( porosity ) * dPorosity * DCO                    
                          * (fvc::grad(CO)&fvc::grad(COa))*CRef
                    + 1.5 * (rmin - rmax)
                          * sqrt( porosity ) * dPorosity * DCR                    
                          * (fvc::grad(CR)&fvc::grad(CRa))*CRef;             
                         
        sensVel =   (U&Ua) * nu * (rmin - rmax) * dInvPermeability; 
               
        
        sens = sensPhi + sensConc + sensVel*correctDim2;


        fvScalarMatrix sensRawEqn
        (
          - fvm::laplacian(Foam::pow(rFilter,2), sensRaw)
          + fvm::Sp(1.0, sensRaw) 
          - sens
          * betaTanh
          //* ((max_gamma-min_gamma)/(brGammaMax - brGammaMin))
          * ( 1.0 - Foam::pow( 
                        Foam::tanh(
                            betaTanh*(   (gammaAux-min_gamma)
                                       / (max_gamma-min_gamma) - etaTanh)
                                     )
                        , 2) 
            )
          / (Foam::tanh(betaTanh*etaTanh) + Foam::tanh(betaTanh*(1.-etaTanh)))
        );
        sensRawEqn.solve();
        
    sensRaw.correctBoundaryConditions();

    forAll(cells, cellI)
    {
        sensRaw[cellI]=meshVol[cellI]*sensRaw[cellI];
        if (centroids[cellI][1] < 0.0)
        {
            sensRaw[cellI]=0.0;
        }
    }


    forAll(cells, cellI)
    {
        // std::cout << "cellI: " << cellI << " sensRaw: " << sensRaw[cellI] << "\n";
        transferVars.brSensitivities[cellI] = sensRaw[cellI];
    }
    // how to copy sensRaw into brSesn...
    //sumSens=gSum(sensRaw);
    //transferVars.brSensitivities[0] = sumSens;
    //transferVars.brSensitivities = sensRaw.cdata();
    // ... copy sensRaw into brSensitivities
    //std::copy(sensRaw,
    //          sensRaw+transferVars.brGammaSize,
    //          transferVars.brSensitivities);
    
    // scalarField as double
    //const double* sensDouble = sensRaw.cdata();
    //transferVars.brSensitivities = sensRaw.cdata();
    //std::copy(sensDouble,
    //          sensDouble + transferVars.brGammaSize,
    //          transferVars.brSensitivities);
}


void newDesignParameters()
{
    
    volScalarField          & rFiberCell        = *transferVars.rFiberCell;
    volScalarField          & invPermeability   = *transferVars.invPermeability;
    volScalarField          & porosity          = *transferVars.porosity;
    volScalarField          & areaPerVol        = *transferVars.areaPerVol;    
    volScalarField          & dInvPermeability  = *transferVars.dInvPermeability;
    volScalarField          & dPorosity         = *transferVars.dPorosity;
    volScalarField          & dAreaPerVol       = *transferVars.dAreaPerVol;
    volScalarField          & gamma             = *transferVars.gamma;
    volScalarField          & gammaAux          = *transferVars.gammaAux;
    volScalarField          & gammaRaw          = *transferVars.gammaRaw;        
    const cellList          & cells             = *transferVars.cells;  
    const dimensionedScalar & rmin              = *transferVars.rmin;
    const dimensionedScalar & rmax              = *transferVars.rmax;
    dimensionedScalar       & unitCellLen       = *transferVars.unitCellLen;
    List<double>            & radiusInterpIn    = *transferVars.radiusInterpIn;
    List<double>            & permInterpIn      = *transferVars.permInterpIn;    
    List<double>            & porosityInterpIn  = *transferVars.porosityInterpIn;
    List<double>            & areaInterpIn      = *transferVars.areaInterpIn;
    List<double>            & DpermInterpIn     = *transferVars.DpermInterpIn;    
    List<double>            & DporosityInterpIn = *transferVars.DporosityInterpIn;
    List<double>            & DareaInterpIn     = *transferVars.DareaInterpIn;       
    const double            & deltaRInterp      = *transferVars.deltaRInterp;
    const volVectorField    & centroids         = *transferVars.centroids;      
    volScalarField          & permeability      = *transferVars.permeability;
    volScalarField          & effDCO            = *transferVars.effDCO;
    volScalarField          & effDCR            = *transferVars.effDCR;
    volScalarField          & effKappa1         = *transferVars.effKappa1;
    volScalarField          & effKappa2         = *transferVars.effKappa2;
    dimensionedScalar       & permeabilityRef   = *transferVars.permeabilityRef;
    dimensionedScalar       & kappa1            = *transferVars.kappa1;
    dimensionedScalar       & kappa2            = *transferVars.kappa2;
    dimensionedScalar       & DCO               = *transferVars.DCO;
    dimensionedScalar       & DCR               = *transferVars.DCR;
    volVectorField          & U                 = *transferVars.U;
    surfaceScalarField      & phi               = *transferVars.phi;
    volScalarField          & p                 = *transferVars.p;
    autoPtr<incompressible::turbulenceModel> & turbulence = *transferVars.turbulence;
    label                   & pRefCell          = *transferVars.pRefCell;
    scalar                  & pRefValue         = *transferVars.pRefValue;
    dimensionedScalar       & nu                = *transferVars.nu;
    fv::options             & fvOptions         = *transferVars.fvOptions;
    volVectorField          & Ua                = *transferVars.Ua;
    surfaceScalarField      & phia              = *transferVars.phia;
    volScalarField          & pa                = *transferVars.pa;
    singlePhaseTransportModel &laminarTransport = *transferVars.laminarTransport; 
    volTensorField          & gradV             = *transferVars.gradV;
    label                   & paRefCell         = *transferVars.paRefCell;        
    scalar                  & paRefValue        = *transferVars.paRefValue;
    double                  & betaTanh          = *transferVars.betaTanh;      
    dimensionedScalar       & etaTanh           = *transferVars.etaTanh;
    dimensionedScalar       & min_gamma         = *transferVars.min_gamma;                                 
    dimensionedScalar       & max_gamma         = *transferVars.max_gamma;
    dimensionedScalar       & rFilter           = *transferVars.rFilter;
    volScalarField          & sensRaw           = *transferVars.sensRaw;
    volScalarField          & bR                = *transferVars.bR;
    volScalarField          & bO                = *transferVars.bO;
    volScalarField          & gR                = *transferVars.gR;
    volScalarField          & gO                = *transferVars.gO;   
    const dimensionedScalar & c1_km             = *transferVars.c1_km;
    const dimensionedScalar & c2_km             = *transferVars.c2_km;
    const dimensionedScalar & c3_km             = *transferVars.c3_km;
    const dimensionedScalar & c4_km             = *transferVars.c4_km;
    const dimensionedScalar & velRef            = *transferVars.velRef;
    const dimensionedScalar & lenRef            = *transferVars.lenRef;     
    volScalarField          & CO                = *transferVars.CO;
    volScalarField          & CR                = *transferVars.CR;
    volScalarField          & Phi1              = *transferVars.Phi1;
    volScalarField          & Phi2              = *transferVars.Phi2;
    volScalarField          & COa               = *transferVars.COa;
    volScalarField          & CRa               = *transferVars.CRa;
    volScalarField          & Psi1              = *transferVars.Psi1;
    volScalarField          & Psi2              = *transferVars.Psi2;
    volScalarField          & kmR               = *transferVars.kmR;
    volScalarField          & kmO               = *transferVars.kmO;    
    dimensionedScalar       & alphaA            = *transferVars.alphaA;
    dimensionedScalar       & alphaC            = *transferVars.alphaC;
    dimensionedScalar       & CRef              = *transferVars.CRef;
    dimensionedScalar       & i0                = *transferVars.i0;
    dimensionedScalar       & FRT               = *transferVars.FRT;
    dimensionedScalar       & nelec             = *transferVars.nelec;
    dimensionedScalar       & stdPotential      = *transferVars.stdPotential;          
    dimensionedScalar       & faradayConstant   = *transferVars.faradayConstant;
    volScalarField          & deltaPhi          = *transferVars.deltaPhi;
    volVectorField          & I1                = *transferVars.I1;
    volVectorField          & I2                = *transferVars.I2;
    dimensionedScalar       & correctDim1UaEqn  = *transferVars.correctDim1UaEqn;
    dimensionedScalar       & correctDim2UaEqn  = *transferVars.correctDim2UaEqn;
    volVectorField          & velDerivInTerm    = *transferVars.velDerivInTerm;
    volVectorField          & cOaGradCO         = *transferVars.cOaGradCO;
    volVectorField          & cRaGradCR         = *transferVars.cRaGradCR;

    /// First, update design variables with output of optimizer
    forAll(cells, cellI)
    {
        gammaRaw[cellI] = transferVars.brGamma[cellI];
        if (centroids[cellI][1] < 0.0)
        {
            gammaRaw[cellI] = 1.0;
        }        
    }
    
    gammaRaw.correctBoundaryConditions();
 
    // Regularize the updated design variable (gammaRaw)     
    fvScalarMatrix gammaAuxEqn
    (
      - fvm::laplacian(Foam::pow(rFilter,2), gammaAux)
      + fvm::Sp(1.0, gammaAux) - gammaRaw              
    );
    gammaAuxEqn.solve();        

    gammaAux = min(max(gammaAux, min_gamma),max_gamma);
    gamma = (max_gamma - min_gamma) *
            (   Foam::tanh(betaTanh*etaTanh)
              + Foam::tanh(betaTanh * (
                                          (gammaAux - min_gamma )
                                        / ( max_gamma - min_gamma )
                                        - etaTanh
                                      )
                          )
            )
            / (   Foam::tanh(betaTanh*etaTanh) 
                + Foam::tanh(betaTanh*(1.-etaTanh))
              )
            + min_gamma;                              
    forAll(cells, cellI)
    {
        if (centroids[cellI][1] < 0.0)
        {
            gamma[cellI]= max_gamma.value();
            gammaRaw[cellI] = 1.0;
        }
    }
    
    gamma.correctBoundaryConditions();
    ///end update of gamma

    
    forAll(cells, cellI)
    {
        double rIn = rmax.value() 
                             + gamma[cellI] 
                             * (rmin.value() - rmax.value());
                               
            int nInterp = floor( ( rIn / unitCellLen.value() 
                                   -radiusInterpIn[0] )
                                 /deltaRInterp );
            
            double r0 = radiusInterpIn[nInterp];
            double r1 = radiusInterpIn[nInterp+1];
            
            double t = ( rIn / unitCellLen.value() - r0 )/
                       ( r1 - r0 );                       

            rFiberCell[cellI] = rIn;     
            invPermeability[cellI] = (  t * permInterpIn[nInterp+1] 
                                     + (1. - t) * permInterpIn[nInterp])
                                     / unitCellLen.value() / unitCellLen.value(); 
            porosity[cellI] = (  t * porosityInterpIn[nInterp+1] 
                                     + (1. - t) * porosityInterpIn[nInterp]);
            areaPerVol[cellI] = (  t * areaInterpIn[nInterp+1] 
                                     + (1. - t) * areaInterpIn[nInterp])
                                     / unitCellLen.value();
            dInvPermeability[cellI] = (  t * DpermInterpIn[nInterp+1] 
                                     + (1. - t) * DpermInterpIn[nInterp])
                                     / unitCellLen.value() / unitCellLen.value()
                                     / unitCellLen.value(); 
            dPorosity[cellI] = (  t * DporosityInterpIn[nInterp+1] 
                                     + (1. - t) * DporosityInterpIn[nInterp])
                                     / unitCellLen.value();
            dAreaPerVol[cellI] = (  t * DareaInterpIn[nInterp+1] 
                                     + (1. - t) * DareaInterpIn[nInterp])
                                     / unitCellLen.value()/unitCellLen.value();
            
            //porosity[cellI] = porosity0.value();
            
            
            /*
            if (gamma[cellI] == max_gamma.value())
            {
                porosity[cellI] = 1.0 - 1.0e-4;
                areaPerVol[cellI] = 0.0;
                invPermeability[cellI] = 0.0;
            }
            */
            
    
            if (centroids[cellI][1] < 0.0)
            {
                porosity[cellI] = 1.0 - 1.0e-4;
                areaPerVol[cellI] = 0.0;
                invPermeability[cellI] = 0.0;
            }
            
            

            if (t > 1.0 && t < 1.000001) {t=1.0;}   

            if ( t  < 0.0 || t > 1.0 )
            {
                Info << "INTERPOLATION ERROR! At rIn: " << rIn << endl;
                Info << "rmin: " << rmin << endl;
                Info << "rmax: " << rmax << endl;
                Info << "radiusInterpIn[0]: " << radiusInterpIn[0] << endl;
                Info << "deltaRInterp: " << deltaRInterp << endl;
                Info << "nInterp: " << nInterp << endl;
                Info << "r0 = " << r0 << endl;
                Info << "r1 = " << r1 << endl;
                Info << "unitCellLen = " << unitCellLen.value() << endl;
            }
        

    }
    
    porosity.correctBoundaryConditions();
    
    

        
        permeability = 1.0 / (invPermeability + 1.0e-20/permeabilityRef);        
        effDCO=sqrt(Foam::pow(porosity,3))*DCO;
        effDCR=sqrt(Foam::pow(porosity,3))*DCR;
        effKappa1=sqrt(Foam::pow(1-porosity,3))*kappa1;
        effKappa2=sqrt(Foam::pow(porosity,3))*kappa2;   
        
        bR = c3_km + c4_km * porosity;
        bO = bR;
            
        gR = DCR/(rFiberCell+1.0e-20*lenRef) 
                                          * c1_km/Foam::pow(porosity+1.0e-20,c2_km)
                                          * Foam::pow(velRef*rFiberCell/DCR,bR);
        gO = DCO/(rFiberCell+1.0e-20*lenRef) 
                                          * c1_km/Foam::pow(porosity+1.0e-20,c2_km)
                                          * Foam::pow(velRef*rFiberCell/DCO,bO);
        gR.correctBoundaryConditions();
        gO.correctBoundaryConditions();
                
        int iterVel = 0;
        double flowConsErr = 1.0;        
                                
        // Pressure-velocity SIMPLE corrector

       
        
        double flowRes = 1.1;
	
        while ( (iterVel < 2000 && flowRes > 1.0e-8) || iterVel < 5)
        // while (iterVel < 2000 && flowRes > 1.0e-8)
        {
            flowRes = 0.0;
            // Momentum predictor
          
            tmp<fvVectorMatrix> tUEqn
            (
                fvm::div(phi, U)
              + turbulence->divDevReff(U)
              + fvm::Sp(nu*invPermeability, U)
             ==
                fvOptions(U)
            );

            fvVectorMatrix& UEqn = tUEqn.ref();



            UEqn.relax();

            fvOptions.constrain(UEqn);
         
            vector spU = solve(UEqn == -fvc::grad(p)).initialResidual();  
            flowRes += spU[0] + spU[1] + spU[2];                               
            
            fvOptions.correct(U);          
            
            volScalarField rAU(1.0/UEqn.A());
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
            tUEqn.clear();
            surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
            adjustPhi(phiHbyA, U, p);

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(p, U, phiHbyA, rAU);

            //You deleted the non-orthogonal corector
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);
                solverPerformance spPress = pEqn.solve();
                
                flowRes += spPress.initialResidual();

                //Info << "Press Res: " << spPress.initialResidual() << endl;

                //if (simple.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }


            // Explicitly relax pressure for momentum corrector
            p.relax();

            // Momentum corrector
            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        

        
            U.storePrevIter();
            p.storePrevIter();
        
                
            iterVel +=1;
            

            
            Info << "Flow Iter: " << iterVel
                 << ", Flow Res: " << flowRes << endl;                                    
        }

        
        Info << "Flow Loop Iter: " << iterVel
                 << ", Flow Loop Res: " << flowRes << endl;
        
        
        // Deleted non-orthogonal corrector 
        {
         
            kmR = gR * Foam::pow(mag(U/velRef)+1.e-20,bR)+1.0e-20*velRef;
            kmO = gO * Foam::pow(mag(U/velRef)+1.e-20,bO)+1.0e-20*velRef;           
                       
            kmR.correctBoundaryConditions();
            kmO.correctBoundaryConditions();


            int iterPot = 0;
            double currConsErr = 1.0; 
            double concConsErr = 1.0;   
            
            double maxScalarsRes = 1.0;
            
            while ( (iterPot < 2000 && maxScalarsRes > 1.0e-8) || iterPot < 5 )
            {  
            
            maxScalarsRes = 0.0;
                                                       
            fvScalarMatrix COEqn
            (
                fvm::div(phi, CO)
              - fvm::laplacian(sqrt(Foam::pow(porosity,3))*DCO, CO)
              + fvm::Sp(
                        areaPerVol*Foam::exp(-1.0*alphaC * FRT * deltaPhi )
                        /
                        (   nelec * faradayConstant * CRef / i0
                          + Foam::exp(alphaA * FRT * deltaPhi ) / kmR
                          + Foam::exp(-1.0*alphaC * FRT * deltaPhi ) / kmO
                        )
                        , CO)
              - areaPerVol * Foam::exp(alphaA * FRT * deltaPhi ) * CR
                        /
                        (   nelec * faradayConstant * CRef / i0
                          + Foam::exp(alphaA * FRT * deltaPhi ) / kmR
                          + Foam::exp(-1.0*alphaC * FRT * deltaPhi ) / kmO
                        )        
            );
                                    
            COEqn.relax();
            solverPerformance spCO = COEqn.solve();
            
            fvScalarMatrix CREqn
            (
                fvm::div(phi, CR)
              - fvm::laplacian(sqrt(Foam::pow(porosity,3))*DCR, CR)
              + fvm::Sp(
                        areaPerVol*Foam::exp(alphaA * FRT * deltaPhi )
                        /
                        (   nelec * faradayConstant * CRef / i0
                          + Foam::exp(alphaA * FRT * deltaPhi ) / kmR
                          + Foam::exp(-1.0*alphaC * FRT * deltaPhi ) / kmO
                        )
                        , CR)
              - areaPerVol * Foam::exp(-1.0*alphaC * FRT * deltaPhi ) * CO
                        /
                        (   nelec * faradayConstant * CRef / i0
                          + Foam::exp(alphaA * FRT * deltaPhi ) / kmR
                          + Foam::exp(-1.0*alphaC * FRT * deltaPhi ) / kmO
                        )     
            );
            
            CREqn.relax();
            solverPerformance spCR = CREqn.solve();
                                                   
           
            fvScalarMatrix Phi1Eqn
            (
                -1.0*fvm::laplacian(sqrt(Foam::pow(1-porosity,3))*kappa1, Phi1)
                + fvm::Sp(areaPerVol * i0 *
                            ( FRT *
                                  (
                                     CR * alphaA 
                                   * Foam::exp(alphaA * FRT * deltaPhi)
                                   + CO * alphaC 
                                   * Foam::exp(-1.0*alphaC * FRT * deltaPhi)
                                  )
                              +   
                                  (
                                     i0 / nelec / CRef * FRT / faradayConstant
                                   * Foam::exp((alphaA-alphaC) * FRT * deltaPhi)
                                   * ( alphaA + alphaC )
                                   * ( CR / kmO + CO / kmR )
                                  )
                             )
                          /
                             Foam::pow(
                                  1. + i0 / nelec / faradayConstant / CRef
                                 *
                                 (  
                                    Foam::exp(alphaA * FRT * deltaPhi ) / kmR
                                  + Foam::exp(-1.0*alphaC * FRT * deltaPhi ) 
                                    / kmO
                                 )
                                 , 2)                            
                           , Phi1)                            
                + areaPerVol * i0 *  
                  (
                     CR * Foam::exp(alphaA * FRT * deltaPhi)
                   - CO * Foam::exp(-1.0*alphaC * FRT * deltaPhi)
                  )
                  /
                  ( 1. + i0/nelec/faradayConstant/CRef
                       * (
                             Foam::exp(alphaA * FRT * deltaPhi ) / kmR
                           + Foam::exp(-1.0*alphaC * FRT * deltaPhi) / kmO
                         )
                  )
                  - areaPerVol * i0 * Phi1 *
                            ( FRT *
                                  (
                                     CR * alphaA 
                                   * Foam::exp(alphaA * FRT * deltaPhi)
                                   + CO * alphaC 
                                   * Foam::exp(-1.0*alphaC * FRT * deltaPhi)
                                  )
                              +   
                                  (
                                     i0 / nelec / CRef * FRT / faradayConstant
                                   * Foam::exp((alphaA-alphaC) * FRT * deltaPhi)
                                   * ( alphaA + alphaC )
                                   * ( CR / kmO + CO / kmR )
                                  )
                             )
                          /
                             Foam::pow(
                                  1. + i0 / nelec / faradayConstant / CRef
                                 *
                                 (  
                                    Foam::exp(alphaA * FRT * deltaPhi ) / kmR
                                  + Foam::exp(-1.0*alphaC * FRT * deltaPhi ) 
                                    / kmO
                                 )
                                 , 2)                           
                            
            );
            
            //Phi1Eqn.relax();
            solverPerformance spPhi1 = Phi1Eqn.solve();
            
            fvScalarMatrix Phi2Eqn
            (
                -1.0*fvm::laplacian(sqrt(Foam::pow(porosity,3))*kappa2, Phi2) 
                + fvm::Sp(areaPerVol * i0 *
                            ( FRT *
                                  (
                                     CR * alphaA 
                                   * Foam::exp(alphaA * FRT * deltaPhi)
                                   + CO * alphaC 
                                   * Foam::exp(-1.0*alphaC * FRT * deltaPhi)
                                  )
                              +   
                                  (
                                     i0 / nelec / CRef * FRT / faradayConstant
                                   * Foam::exp((alphaA-alphaC) * FRT * deltaPhi)
                                   * ( alphaA + alphaC )
                                   * ( CR / kmO + CO / kmR )
                                  )
                             )
                          /
                             Foam::pow(
                                  1. + i0 / nelec / faradayConstant / CRef
                                 *
                                 (  
                                    Foam::exp(alphaA * FRT * deltaPhi ) / kmR
                                  + Foam::exp(-1.0*alphaC * FRT * deltaPhi ) 
                                    / kmO
                                 )
                                 , 2)                            
                           , Phi2)                            
                - areaPerVol * i0 *  
                  (
                     CR * Foam::exp(alphaA * FRT * deltaPhi)
                   - CO * Foam::exp(-1.0*alphaC * FRT * deltaPhi)
                  )
                  /
                  ( 1. + i0/nelec/faradayConstant/CRef
                       * (
                             Foam::exp(alphaA * FRT * deltaPhi ) / kmR
                           + Foam::exp(-1.0*alphaC * FRT * deltaPhi) / kmO
                         )
                  )
                  - areaPerVol * i0 * Phi2 *
                            ( FRT *
                                  (
                                     CR * alphaA 
                                   * Foam::exp(alphaA * FRT * deltaPhi)
                                   + CO * alphaC 
                                   * Foam::exp(-1.0*alphaC * FRT * deltaPhi)
                                  )
                              +   
                                  (
                                     i0 / nelec / CRef * FRT / faradayConstant
                                   * Foam::exp((alphaA-alphaC) * FRT * deltaPhi)
                                   * ( alphaA + alphaC )
                                   * ( CR / kmO + CO / kmR )
                                  )
                             )
                          /
                             Foam::pow(
                                  1. + i0 / nelec / faradayConstant / CRef
                                 *
                                 (  
                                    Foam::exp(alphaA * FRT * deltaPhi ) / kmR
                                  + Foam::exp(-1.0*alphaC * FRT * deltaPhi ) 
                                    / kmO
                                 )
                                 , 2)       
            );
            
            //Phi2Eqn.relax();
            solverPerformance spPhi2 = Phi2Eqn.solve();
            
            deltaPhi = Phi1 -Phi2 - stdPotential;

            I2 = -1.0*sqrt(Foam::pow(porosity,3))*kappa2*fvc::grad(Phi2);
            I1 = -1.0*sqrt(Foam::pow(1.-porosity,3))*kappa1*fvc::grad(Phi1);    
            
            
            maxScalarsRes += spCO.initialResidual() + spCR.initialResidual()
                         + spPhi1.initialResidual() + spPhi2.initialResidual();

            iterPot +=1;

            Info << "Iter Scalars: " << iterPot  
                 <<", Sum of Scalars Residuals: " << maxScalarsRes 
                 << endl; 
                                                            
            }
            
            Info << "Scalar Loop Iter: " << iterPot  
                 <<", Scalar Loop Res: " << maxScalarsRes 
                 << endl; 
            
            // Outer adjnts while loop
            int iterAdjnts = 0;
            double maxAdjRes = 1.0; 
            
            while ( (iterAdjnts < 2000 && maxAdjRes > 1.0e-8) || iterAdjnts < 5)
            {            
            
            maxAdjRes = 0.0;
            
            fvScalarMatrix COaEqn
            (
              - fvm::div(phi, COa)
              - fvm::laplacian(sqrt(Foam::pow(porosity,3))*DCO, COa)
              + fvm::Sp(areaPerVol * Foam::exp(-1.0*alphaC * FRT * deltaPhi)
                        /
                        (
                           nelec * faradayConstant * CRef / i0
                         + Foam::exp(alphaA * FRT * deltaPhi ) / kmR
                         + Foam::exp(-1.0*alphaC * FRT * deltaPhi ) / kmO
                        )
                        , COa)              
              - areaPerVol * nelec * faradayConstant
                        * Foam::exp(-1.0*alphaC * FRT * deltaPhi)
                        /
                        (
                           nelec * faradayConstant * CRef / i0
                         + Foam::exp(alphaA * FRT * deltaPhi ) / kmR
                         + Foam::exp(-1.0*alphaC * FRT * deltaPhi ) / kmO
                        )
                        * ( CRa / nelec / faradayConstant + Psi1 - Psi2)
            );
            
            COaEqn.relax();
            solverPerformance spCOa = COaEqn.solve();            
            
            fvScalarMatrix CRaEqn
            (
              - fvm::div(phi, CRa)
              - fvm::laplacian(sqrt(Foam::pow(porosity,3))*DCR, CRa)
              + fvm::Sp(areaPerVol * Foam::exp(alphaA * FRT * deltaPhi)
                        /
                        (
                           nelec * faradayConstant * CRef / i0
                         + Foam::exp(alphaA * FRT * deltaPhi ) / kmR
                         + Foam::exp(-1.0*alphaC * FRT * deltaPhi ) / kmO
                        )
                        , CRa)              
              + areaPerVol * nelec * faradayConstant
                        * Foam::exp(alphaA * FRT * deltaPhi)
                        /
                        (
                           nelec * faradayConstant * CRef / i0
                         + Foam::exp(alphaA * FRT * deltaPhi ) / kmR
                         + Foam::exp(-1.0*alphaC * FRT * deltaPhi ) / kmO
                        )
                        * ( -1.*COa / nelec / faradayConstant + Psi1 - Psi2)
            );
            
            CRaEqn.relax();
            solverPerformance spCRa = CRaEqn.solve();
                                
            
            fvScalarMatrix Psi1Eqn
            (
                -1.0*fvm::laplacian(sqrt(Foam::pow(1.-porosity,3))*kappa1, Psi1) 
                +fvm::Sp(areaPerVol * i0 *
                            ( FRT *
                                  (
                                     CR * alphaA 
                                   * Foam::exp(alphaA * FRT * deltaPhi)
                                   + CO * alphaC 
                                   * Foam::exp(-1.0*alphaC * FRT * deltaPhi)
                                  )
                              +   
                                  (
                                     i0 / nelec / CRef * FRT / faradayConstant
                                   * Foam::exp((alphaA-alphaC) * FRT * deltaPhi)
                                   * ( alphaA + alphaC )
                                   * ( CR / kmO + CO / kmR )
                                  )
                             )
                          /
                             Foam::pow(
                                  1. + i0 / nelec / faradayConstant / CRef
                                 *
                                 (  
                                    Foam::exp(alphaA * FRT * deltaPhi ) / kmR
                                  + Foam::exp(-1.0*alphaC * FRT * deltaPhi ) 
                                    / kmO
                                 )
                                 , 2)                            
                           , Psi1)  
                           
                           
                 + (areaPerVol * i0 *
                            ( FRT *
                                  (
                                     CR * alphaA 
                                   * Foam::exp(alphaA * FRT * deltaPhi)
                                   + CO * alphaC 
                                   * Foam::exp(-1.0*alphaC * FRT * deltaPhi)
                                  )
                              +   
                                  (
                                     i0 / nelec / CRef * FRT / faradayConstant
                                   * Foam::exp((alphaA-alphaC) * FRT * deltaPhi)
                                   * ( alphaA + alphaC )
                                   * ( CR / kmO + CO / kmR )
                                  )
                             )
                          /
                             Foam::pow(
                                  1. + i0 / nelec / faradayConstant / CRef
                                 *
                                 (  
                                    Foam::exp(alphaA * FRT * deltaPhi ) / kmR
                                  + Foam::exp(-1.0*alphaC * FRT * deltaPhi ) 
                                    / kmO
                                 )
                                 , 2)                            
                    )*( (CRa-COa) / nelec / faradayConstant - Psi2)            
                                                                                                                            
            );
            solverPerformance spPsi1 = Psi1Eqn.solve();            
            
            fvScalarMatrix Psi2Eqn
            (
                -1.0*fvm::laplacian(sqrt(Foam::pow(porosity,3))*kappa2, Psi2) 
                +fvm::Sp(areaPerVol * i0 *
                            ( FRT *
                                  (
                                     CR * alphaA 
                                   * Foam::exp(alphaA * FRT * deltaPhi)
                                   + CO * alphaC 
                                   * Foam::exp(-1.0*alphaC * FRT * deltaPhi)
                                  )
                              +   
                                  (
                                     i0 / nelec / CRef * FRT / faradayConstant
                                   * Foam::exp((alphaA-alphaC) * FRT * deltaPhi)
                                   * ( alphaA + alphaC )
                                   * ( CR / kmO + CO / kmR )
                                  )
                             )
                          /
                             Foam::pow(
                                  1. + i0 / nelec / faradayConstant / CRef
                                 *
                                 (  
                                    Foam::exp(alphaA * FRT * deltaPhi ) / kmR
                                  + Foam::exp(-1.0*alphaC * FRT * deltaPhi ) 
                                    / kmO
                                 )
                                 , 2)                            
                           , Psi2)  
                           
                           
                 - (areaPerVol * i0 *
                            ( FRT *
                                  (
                                     CR * alphaA 
                                   * Foam::exp(alphaA * FRT * deltaPhi)
                                   + CO * alphaC 
                                   * Foam::exp(-1.0*alphaC * FRT * deltaPhi)
                                  )
                              +   
                                  (
                                     i0 / nelec / CRef * FRT / faradayConstant
                                   * Foam::exp((alphaA-alphaC) * FRT * deltaPhi)
                                   * ( alphaA + alphaC )
                                   * ( CR / kmO + CO / kmR )
                                  )
                             )
                          /
                             Foam::pow(
                                  1. + i0 / nelec / faradayConstant / CRef
                                 *
                                 (  
                                    Foam::exp(alphaA * FRT * deltaPhi ) / kmR
                                  + Foam::exp(-1.0*alphaC * FRT * deltaPhi ) 
                                    / kmO
                                 )
                                 , 2)                            
                    )*( (CRa-COa) / nelec / faradayConstant + Psi1)
            );               
            solverPerformance spPsi2 = Psi2Eqn.solve(); 
            
            maxAdjRes += spCOa.initialResidual() + spCRa.initialResidual()
                         + spPsi1.initialResidual() + spPsi2.initialResidual();

            iterAdjnts += 1;                         
            Info << "Iter Adnt: " << iterAdjnts << 
            ", Sum of Adjnt Initial Residuals: " << maxAdjRes << endl;                    
            
            } 
            
            Info << "Scalar Adj Loop Iter: " << iterAdjnts << 
            ", Scalar Adj Loop Res: " << maxAdjRes << endl;
            
        }
                 
        // Outer velocity adjnts while loop
        
        Ua.storePrevIter();
        pa.storePrevIter();
         
        
       int iterVelAdjnts = 0;
       double maxVelAdjRes = 1.0; 
        
        while ( (iterVelAdjnts < 2000 && maxVelAdjRes > 1.0e-8) || iterVelAdjnts < 5)     
        {

            maxVelAdjRes = 0.0;            
            
            // Adjoint Momentum predictor


            velDerivInTerm =
                   (                  
                      areaPerVol * nelec * faradayConstant
                    * (
                          CO * Foam::exp(-1.0*alphaC * FRT * deltaPhi) 
                        - CR * Foam::exp(alphaA * FRT * deltaPhi)
                      ) 
                    /
                      Foam::pow(
                          nelec * faradayConstant * CRef / i0
                        + Foam::exp(-1.0*alphaC * FRT * deltaPhi)
                          / (gO * Foam::pow(mag(U/velRef)+1.e-20,bO))
                        + Foam::exp(alphaA * FRT * deltaPhi)
                          / (gR * Foam::pow(mag(U/velRef)+1.e-20,bR))
                         ,2)
                    * (
                          Foam::exp(alphaA * FRT * deltaPhi)*bR
                          / (gR * Foam::pow(mag(U/velRef)+1.e-20,bR+1.))
                        + Foam::exp(-1.0*alphaC * FRT * deltaPhi)*bO
                          / (gO * Foam::pow(mag(U/velRef)+1.e-20,bO+1.))
                      )
                   * ( (CRa-COa) /nelec / faradayConstant + Psi1 - Psi2 )
                   * correctDim2UaEqn
                   * U/(mag(U) + 1.e-100*velRef)                  
                   );

            cOaGradCO = COa*correctDim1UaEqn*fvc::grad(CO);
            cRaGradCR = CRa*correctDim1UaEqn*fvc::grad(CR);    
            gradV = fvc::grad(U);  

            tmp<fvVectorMatrix> tUaEqn
            (
                fvm::div(-phi, Ua)
              - fvm::laplacian(nu, Ua)
              + fvm::Sp(nu*invPermeability, Ua)
              + fvm::SuSp(tr(gradV)/3.0,Ua)
              + (dev(gradV) & Ua)
              + cOaGradCO + cRaGradCR - velDerivInTerm              
            );
            fvVectorMatrix& UaEqn = tUaEqn.ref();

            UaEqn.relax();

            fvOptions.constrain(UaEqn);    


            vector spUa = solve(UaEqn == -fvc::grad(pa)).initialResidual();  
            maxVelAdjRes += spUa[0] + spUa[1] + spUa[2];

            fvOptions.correct(Ua);
            
            volScalarField rAUa(1.0/UaEqn.A());
            volVectorField HbyAa("HbyAa", Ua);
            HbyAa = rAUa*UaEqn.H();
            tUaEqn.clear();
            surfaceScalarField phiHbyAa("phiHbyAa", fvc::flux(HbyAa));
            adjustPhi(phiHbyAa, Ua, pa);

            // Non-orthogonal pressure corrector loop
            {
                fvScalarMatrix paEqn
                (
                    fvm::laplacian(rAUa, pa) == fvc::div(phiHbyAa)
                );

                paEqn.setReference(paRefCell, paRefValue);
                
                solverPerformance spPressAdj = paEqn.solve();
                
                maxVelAdjRes += spPressAdj.initialResidual();

                //if (simple.finalNonOrthogonalIter())
                {
                    phia = phiHbyAa - paEqn.flux();
                }
            }

            // Explicitly relax pressure for adjoint momentum corrector
            pa.relax();

            // Adjoint momentum corrector
            Ua = HbyAa - rAUa*fvc::grad(pa);
            Ua.correctBoundaryConditions();
            //fvOptions.correct(Ua);
            
            Ua.storePrevIter();
            pa.storePrevIter();
        
                
            iterVelAdjnts +=1;
            
            Info << "Flow Adj Iter: " << iterVelAdjnts
                 << ", Adj Flow Res: " << maxVelAdjRes << endl;
                 
        }
        
        Info << "Flow Adj Loop Iter: " << iterVelAdjnts
                 << ", Flow Adj Loop Res: " << maxVelAdjRes << endl;

        laminarTransport.correct();         
    
} // newDesignParameters
// ************************************************************************* //
