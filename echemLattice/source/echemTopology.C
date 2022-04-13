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
#include "PstreamGlobals.H"
#include "IFstream.H"
#include <iomanip>

#include <map>
#include <vector>
#include <algorithm>
#include <functional>
#include <iomanip>

#include "op.hpp"
#include "op_debug.hpp"
#include "nlopt_op.hpp"

// rfunction 
#include "rfunc.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
static const double solver_tol = 1.0e-8;
static const long interation_limit = 5000;
static const double starting_value = 0.2;
static const double constraint_value = 0.3;  // Upper volume fraction limit
// constraint_value turned of if >= 1.0
static const double gamma_raw_min_value = 0.1;  // This is brGammaMin (design variable min)
static const double gamma_raw_max_value = 0.4;  // This is brGammaMax (design variable max)
// make l_design_domain negative to turn off
static const double l_design_domain = -1e-3;  // y_max minus this value is top domain
// number of cells in one dual mesh cell
static const long n_cells_into_dual_mesh = 20.0;

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
    dimensionedScalar * bR;
    dimensionedScalar * bO;
    dimensionedScalar * gR;
    dimensionedScalar * gO;
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
    Foam::fvMesh * mesh;
    std::unordered_map<long, bool> * cell_can_change;  // cell_can_change[cellI] is True when the cell is allowed to change
    std::unordered_map<long, bool> * is_top_cell;  //is_top_cell[cellI] is True when it's a top cell!
    int * number_of_design_var;  // a new design variable counter that doesn't mess up the older ones
    std::unordered_map<long, std::vector<long>> * global_var_to_local_cell_id;
    std::unordered_map<long, long> * local_cell_id_to_global_var;
    volScalarField * gammaId;
    volScalarField * cellChangeMask;
    volScalarField * gammaIdx;
    volScalarField * gammaIdz;
    volScalarField * rank;
    volScalarField * rank_label;
    double * qArea;
    double * qPerm; 
    dimensionedScalar * a0;
    dimensionedScalar * porosity0;
    dimensionedScalar * permeability0;    
    volScalarField * gammaGeom;
    volScalarField * gammaGeomGrad;
    volScalarField * gammaGeomPlus;
    volScalarField * gammaGeomMinus;
}; 

struct dataTransfer transferVars;


void setNumberDecisionVars();

void setStartingValue();        

void setBounds();        

void saveResults();

void saveResults(std::function<void()> &);

double evalObjective();

double evalConstraint();

void evalConstraintSensitivity();

void newDesignParameters();

void evalSensitivities();

void setCellCanChange(auto & cell_can_change)
{
    const cellList          & cells             = *transferVars.cells;
    const volVectorField    & centroids         = *transferVars.centroids; 

    forAll(cells, cellI)
    {
      cell_can_change[cellI] = true;
      if (centroids[cellI][1] < 0.0)
      {
        cell_can_change[cellI] = false;
      }
    }
}

void setCellCanChange2dDomain(double global_y_max)
{
    // this should be called after the global dimensions of the bounding box
    const cellList          & cells             = *transferVars.cells;
    Foam::fvMesh            & mesh              = *transferVars.mesh;
    auto & cell_can_change = *transferVars.cell_can_change;
    auto & is_top_cell = *transferVars.is_top_cell;
    auto & cellChangeMask = *transferVars.cellChangeMask;

    forAll(cells, cellI)
    {
      is_top_cell[cellI] = false;
      if (cell_can_change[cellI] )
      {
        auto Y_pos_local = mesh.C()[cellI].y();
        if (global_y_max - l_design_domain <= Y_pos_local)
        {
          is_top_cell[cellI] = true;
          cell_can_change[cellI] = false;
        }
      }
	    cellChangeMask[cellI] = cell_can_change[cellI];
    }
}

std::tuple<double, double, double> getCellDimensions() {
    const cellList          & cells             = *transferVars.cells;
    Foam::fvMesh            & mesh              = *transferVars.mesh;
    auto     & cell_can_change      = *transferVars.cell_can_change;
    const faceList & ff = mesh.faces();
    const pointField & pp = mesh.points();

    // Need to loop through all cells and find the first changeable cell
    long cell_index = 0;
    bool a_cell_can_change = false;
    forAll(cells, cellI)
    {
      if (cell_can_change[cellI])
      {
        cell_index = cellI;
        a_cell_can_change = true;
        break;
      }
    }
    const cell & cc = mesh.cells()[cell_index];
    labelList pLabels(cc.labels(ff));
    pointField pLocal(pLabels.size(), vector::zero);

    forAll (pLabels, pointi)
    {
	    pLocal[pointi] = pp[pLabels[pointi]];
    }
    if (a_cell_can_change){
      return std::make_tuple(
        Foam::max(pLocal & vector(1,0,0)) - Foam::min(pLocal & vector(1,0,0)), // xdim
        Foam::max(pLocal & vector(0,1,0)) - Foam::min(pLocal & vector(0,1,0)), // ydim
        Foam::max(pLocal & vector(0,0,1)) - Foam::min(pLocal & vector(0,0,1)) // zdim
        );
    }
    else {
      double xdim = HUGE_VAL;
      double ydim = HUGE_VAL;
      double zdim = HUGE_VAL;
      return std::make_tuple(xdim, ydim, zdim); 
    }

}

std::tuple<double, double, double, double, double, double> getBoundingBox() {
    double x_min = HUGE_VAL;
    double y_min = HUGE_VAL;
    double z_min = HUGE_VAL;
    double x_max = -HUGE_VAL;
    double y_max = -HUGE_VAL;
    double z_max = -HUGE_VAL;
    const cellList          & cells             = *transferVars.cells;
    auto     & cell_can_change      = *transferVars.cell_can_change;
    Foam::fvMesh            & mesh              = *transferVars.mesh;
  
    forAll(cells, cellI)
    {
      if (cell_can_change[cellI])
      {
        auto X_pos_local = mesh.C()[cellI].x();
        auto Y_pos_local = mesh.C()[cellI].y();
        auto Z_pos_local = mesh.C()[cellI].z();
        x_min = min(X_pos_local, x_min);
        x_max = max(X_pos_local, x_max);
        y_min = min(Y_pos_local, y_min);
        y_max = max(Y_pos_local, y_max);
        z_min = min(Z_pos_local, z_min);
        z_max = max(Z_pos_local, z_max);
      }
    }
    return std::make_tuple(x_min, y_min, z_min, x_max, y_max, z_max);
}

int main(int argc, char *argv[])
{

  using ADType = autodiff::real;
  
    Foam::argList args(argc, argv); 

    
    if (!args.checkRootCase())
    {
        Foam::FatalError.exit();
    } 

    
    
    Foam::Info<< "Create time\n" << Foam::endl;
    Foam::Time runTime(Foam::Time::controlDictName, args);

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
                                 dimless, 1e-2);
    dimensionedScalar max_gamma("max_gamma",
                                 dimless, 0.99);

//    dimensionedScalar min_gammaRaw("min_gammaRaw",
//                                 dimless, 0.0);
//    dimensionedScalar max_gammaRaw("max_gammaRaw",
//                                 dimless, 1.0);

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
    //double qArea = qArea0_.value();    
    //double qPerm = qPermInterp0_.value();
    double gammaTrnPerm = gammaTrnPermInterp0_.value();   

    Info << "qArea: " << qArea0_.value() << endl;        
    Info << "qPerm: " << qPermInterp0_.value() << endl;
    Info << "betaTanh: " << betaTanh << endl;  
    Info << "gammaTrnPerm: " << gammaTrnPerm << endl;      
    
    // Variale step-length parameter
    
    double maxStepLength = stepLength.value();
    int restepFlag = 0;
    
    stepLength = stepLengthInit;


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
    double VT = constraint_value; 
    int tempNumbDesignVars = brGammaSize;
    double dv_minBr=gamma_raw_min_value;
    double dv_maxBr=gamma_raw_max_value;


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
        new double [brGammaSize], //brGammaStart, empty array init
        &brGammaSize,
        &dv_minBr, //brGammaMin, global min value of design variable 
        &dv_maxBr, //brGammaMax, global max value of design variable        
        new double [brGammaSize], //brGamma, empty array init
        new double [brGammaSize], //brSensitivities, empty array init
        new double [brGammaSize], //brLowerBound, empty array init
        new double [brGammaSize], //brUpperBound, empty array init
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
        &VT,
        &mesh,  // added mesh
        new std::unordered_map<long, bool>(), // cell_can_change
        new std::unordered_map<long, bool>(), // is_top_cell
        &tempNumbDesignVars, // number_of_design_var
        new std::unordered_map<long, std::vector<long>>(), // global_var_to_local_cell_id
        new std::unordered_map<long, long>(), //local_cell_id_to_global_var;
        &gammaId,
        &cellChangeMask,
        &gammaIdx,
        &gammaIdz,
        &rank,
        &rank_label,
        &qArea0_.value(),
        &qPermInterp0_.value(),
        &a0,
        &porosity0,
        &permeability0,
        &gammaGeom,
        &gammaGeomGrad,
        &gammaGeomPlus,
        &gammaGeomMinus
    };

    U.storePrevIter();
    p.storePrevIter();
   
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

    int my_rank = op::mpi::getRank(mpicomm);
    int size = op::mpi::getNRanks(mpicomm);  // total number of procs
    bool no_design_vars = true;

    // This populates the maps which are used to determine which OF cell ID belongs to which
    // design variable.

    long x_cells, y_cells, z_cells;
    double global_x_min, global_y_min, global_z_min;
    { // Formerly setBounds()
      const cellList          & cells             = *transferVars.cells;
      auto & cell_can_change = *transferVars.cell_can_change;
      auto & global_var_to_local_cell_id = *transferVars.global_var_to_local_cell_id;
      auto & local_cell_id_to_global_var = *transferVars.local_cell_id_to_global_var;
      // find the inital cells that can change
      setCellCanChange(cell_can_change);

      // get the bounding box for this partition
      auto [x_min, y_min, z_min, x_max, y_max, z_max] = getBoundingBox();

      // get global bounding box
      double global_x_max, global_y_max, global_z_max;

      op::mpi::Allreduce(x_min, global_x_min, MPI_MIN, mpicomm);
      op::mpi::Allreduce(y_min, global_y_min, MPI_MIN, mpicomm);
      op::mpi::Allreduce(z_min, global_z_min, MPI_MIN, mpicomm);

      op::mpi::Allreduce(x_max, global_x_max, MPI_MAX, mpicomm);
      op::mpi::Allreduce(y_max, global_y_max, MPI_MAX, mpicomm);
      op::mpi::Allreduce(z_max, global_z_max, MPI_MAX, mpicomm);

      // update cells based on their Y position of designable space
      setCellCanChange2dDomain(global_y_max);  // needs to be before xDim yDim zDim calc!

      // get changeable cell dimensions
      auto [localxDim, localyDim, localzDim] = getCellDimensions();
      double xDim;
      double yDim;
      double zDim;

      // this is a work arround to a bug that happens when there are no changeable cells on a rank
      // but essentially one of these dimmensions ends up too big
      op::mpi::Allreduce(localxDim, xDim, MPI_MIN, mpicomm);
      op::mpi::Allreduce(localyDim, yDim, MPI_MIN, mpicomm);
      op::mpi::Allreduce(localzDim, zDim, MPI_MIN, mpicomm);

      // int total_cells = cells.size();
      double len_x = global_x_max - global_x_min;
      double len_y = global_y_max - global_y_min;
      double len_z = global_z_max - global_z_min;
      // need to add 1 because looking at centers
      std::cout << "  : len_x: " << len_x << " len_y: " << len_y << " len_z: " << len_z << "\n";
      std::cout << "  : xDim: " << xDim << " yDim: " << yDim << " zDim: " << zDim << " MyRank: " << my_rank << "\n";

      // center edge to center edge misses 1 element
      x_cells = std::lround((len_x / xDim + 0.5) / n_cells_into_dual_mesh);
      y_cells = std::lround((len_y / yDim + 0.5) / n_cells_into_dual_mesh);
      z_cells = std::lround((len_z / zDim + 0.5) / n_cells_into_dual_mesh);
      // These must be the same on each rank
      std::cout << "  : x_cells: " << x_cells << " y_cells: " << y_cells << " z_cells: " << z_cells << "\n";
      // design variables are columns in Y, and the plane is the XY plane
      forAll(cells, cellI)
      {
        if (cell_can_change[cellI])
        {
          double x = mesh.C()[cellI].x();
          double y = mesh.C()[cellI].y();
          double z = mesh.C()[cellI].z();

          long x_int = std::floor((x-global_x_min) / xDim + 0.5);
          long y_int = std::floor((y-global_y_min) / yDim + 0.5);
          long z_int = std::floor((z-global_z_min) / zDim + 0.5);

          x_int = std::floor(x_int/n_cells_into_dual_mesh + (1/n_cells_into_dual_mesh) - (1/(n_cells_into_dual_mesh*2)));
          y_int = std::floor(y_int/n_cells_into_dual_mesh + (1/n_cells_into_dual_mesh) - (1/(n_cells_into_dual_mesh*2)));
          z_int = std::floor(z_int/n_cells_into_dual_mesh + (1/n_cells_into_dual_mesh) - (1/(n_cells_into_dual_mesh*2)));

          // design variable have a unique x & z value     
	  std::size_t label = (x_int*z_cells*y_cells) + (z_int*y_cells) + y_int; 
          global_var_to_local_cell_id[label].push_back(cellI);
          local_cell_id_to_global_var[cellI] = label;
          gammaId[cellI] = label;
          gammaIdx[cellI] = x_int;
          gammaIdz[cellI] = z_int;
          no_design_vars = false;
        }
	rank[cellI] = my_rank;
      }
    }


    
    auto & global_var_to_local_cell_id = *transferVars.global_var_to_local_cell_id;
    std::vector<std::size_t> global_ids_on_rank;  // this is the keys for global_var_to_local_cell_id

    for (auto & [key, _] : global_var_to_local_cell_id) {
      global_ids_on_rank.push_back(key);
    }
    

    std::cout << "rank initial: " << my_rank << " ";
    std::unordered_map<std::size_t, std::size_t> label_index_map;  // map a label to an index vector
    std::unordered_map<std::size_t, std::size_t> index_label_map;  // map an index vector to a label
    
    // design variables on a rank are ordered from 0, 1, ..., to ndvar
    std::stringstream label_file_name;
    label_file_name << "pattern." << my_rank;
    std::ofstream label_file(label_file_name.str());
    {
      int id_counter = 0;
      for (auto id : global_ids_on_rank) {
        label_file << id << " ";
        label_index_map[id] = id_counter;
        index_label_map[id_counter] = id;
        id_counter++;
      }
      label_file << std::endl;
    }
    
    

    // print out label corresponding to cell id (debugging)
    {
      auto & cell_can_change = *transferVars.cell_can_change;
      auto & local_cell_id_to_global_var = *transferVars.local_cell_id_to_global_var;
      const cellList          & cells             = *transferVars.cells;
      forAll(cells, cellI)
      {
        if (cell_can_change[cellI])
        {
          auto lbl = local_cell_id_to_global_var[cellI];
          auto idx = label_index_map[lbl];
          rank_label[cellI] = global_ids_on_rank[idx];
        }
      }
    }
    constexpr int root = 0;

    auto comm_pattern= op::AdvancedRegistration(global_ids_on_rank, root, mpicomm);

    // For debugging print out comm_pattern
    ////op::debug::writeCommPatternToDisk(comm_pattern, my_rank);
    /** End Registration **/
    
    std::fill_n (transferVars.brSensitivities, brGammaSize, 0.0);

    int & number_of_design_var = *transferVars.number_of_design_var;
    int actual_y_cols = global_ids_on_rank.size();
    std::cout << "  : actual y cols: " << actual_y_cols << " my_rank: " << my_rank << "\n";
    std::vector<double> x(actual_y_cols); //Defines local design variables


    for (auto & dv : x) {
      dv = starting_value;  // const defined below includes
    }
    setStartingValue();  // Need to copy/paste code from update into setStartingValue for geom

    std::cout << "  : setStartingValue Done! rank: " << my_rank << "\n";

    // Define op primitives
    using opVectorType = std::vector<double>;
    auto lowerBoundFn = [&]() { return std::vector<double>(actual_y_cols, *transferVars.brGammaMin);};
    auto upperBoundFn = [&]() { return std::vector<double>(actual_y_cols, *transferVars.brGammaMax);};
    op::Vector<opVectorType> op_variables(x,  lowerBoundFn, upperBoundFn);

    auto nlopt_options = op::NLoptOptions{.Int = {{"maxeval", 10000}},
                                          .Double = {{"xtol_rel", 1e-16},
                                                     {"verbosity", 2}/*,
                                                     {"inner_maxeval", 10}*/},
                                          .String = {{}},
                                          .algorithm = nlopt::LD_MMA};  // choices are LD_MMA, LD_CCSAQ, & LD_SLSQP

    auto opt = op::NLopt(op_variables, nlopt_options, mpicomm, comm_pattern);

    double scale = 1.0; // objective scaling shouldn't really be needed on the echem problem

    // Define local objective function
    auto global_objective = [&]([[maybe_unused]] const opVectorType & local_variables) {
      double obj = evalObjective();
      return obj * scale;
    };
  
    // Local integration/summing of objective sensitivities
    auto local_obj_grad = [&]([[maybe_unused]] const opVectorType & local_variables) {
      opVectorType grad(local_variables.size(), 0.); // initialize gradients as zero
      evalSensitivities();  // eval sensitivities...
      // copied from EvaluateObjectivetGradient
      auto & global_var_to_local_cell_id = *transferVars.global_var_to_local_cell_id;

      // loop through all columns..
      for (size_t i = 0; i < grad.size(); i++)
      {
        size_t global_design_var = global_ids_on_rank[i];  // this is the index of the 'global design variable'
        // std::cout << "\n  : global_des_var: " << global_design_var << " Index label map: " << index_label_map[i] << " \n";
        auto n_cells = global_var_to_local_cell_id[global_design_var].size();
        for (size_t j = 0; j < n_cells; j++)
        {
          // add the local cell contribution to the sens
          size_t local_cell = global_var_to_local_cell_id[global_design_var][j];
          grad[i] += transferVars.brSensitivities[local_cell];
        }
        grad[i] *= scale;
        //std::cout << '  : grad[' << i << '] = ' << grad[i] << ' ,  x[' << i << '] = ' << x[i] << std::endl;
        std::cout <<  grad[i] << std::endl;
      }
      
      for (size_t i = 0; i < grad.size(); i++)
      {
        std::cout <<  x[i] << std::endl;
      }

      return grad;
    };

    auto reduced_local_obj_grad =
      opt.generateReducedLocalGradientFunction(local_obj_grad,
                op::utility::reductions::sumOfCollection<std::vector<double>>);
    
    // op needs reduced_local_obj_gradient
    op::Functional obj(global_objective, reduced_local_obj_grad);
    
    auto global_constraint = [&]([[maybe_unused]] const opVectorType & local_variables) {
      double constraint = evalConstraint();    
      return constraint;
    };

    // Local integration/summing of constraint sensitivities
    auto local_constraint_grad = [&]([[maybe_unused]] const opVectorType & local_variables) {
      opVectorType grad(local_variables.size(), 0.); // initialize gradients as zero
      // copied from EvaluateConstraintGradients
      evalConstraintSensitivity();
      auto & global_var_to_local_cell_id = *transferVars.global_var_to_local_cell_id;
      // We only have one constraint
      // loop through all columns..
      for (size_t i = 0; i < grad.size(); i++)
      {
        // auto n_cells = global_var_to_local_cell_id[i].size();
        auto n_cells = global_var_to_local_cell_id[global_ids_on_rank[i]].size();
        for (size_t j = 0; j < n_cells; j++)
        {
          // add the local cell contribution to the sens
          // grad[i] += transferVars.VCsens[global_var_to_local_cell_id[i][j]];
          grad[i] += transferVars.VCsens[global_var_to_local_cell_id[global_ids_on_rank[i]][j]];

        }
      }
      return grad;
    };

    auto reduced_local_constraint_grad =
      opt.generateReducedLocalGradientFunction(local_constraint_grad,
                op::utility::reductions::sumOfCollection<std::vector<double>>);

    double cons1_lb = 0.; // lower box constraint
    double cons1_ub = *transferVars.VT;  //upper box constraint
    op::Functional cons1(global_constraint, reduced_local_constraint_grad, cons1_lb, cons1_ub);

    
    //auto unit_cell = SmoothHeaviside(k, sphere);
    
    ////////////////////////////////////////////////////////////////////////////
    
    // what to do when we get new decision variables
    

    double nCell_rod = 3.0;
    double k_heaviside = 1.0;

    opt.update = [&, actual_y_cols, x_cells, y_cells, z_cells]() { // copied from NewXCallback()

      // We need to update transferVars.brGamma given these new dv
      auto    & local_cell_id_to_global_var = *transferVars.local_cell_id_to_global_var;
      auto     & cell_can_change      = *transferVars.cell_can_change;
      Foam::fvMesh            & mesh              = *transferVars.mesh;
      volScalarField          & gammaGeom             = *transferVars.gammaGeom;
      volScalarField          & gammaGeomGrad             = *transferVars.gammaGeomGrad;
      volScalarField          & gammaGeomPlus             = *transferVars.gammaGeomPlus;
      volScalarField          & gammaGeomMinus             = *transferVars.gammaGeomMinus;
      
      // Compose lattice function
    Point3D<ADType> pmx; // point minus x
    Point3D<ADType> ppx; // point plus  x
    Point3D<ADType> pmy; // point minus y
    Point3D<ADType> ppy; // point plus  y
    Point3D<ADType> pmz; // point minus z
    Point3D<ADType> ppz; // point plus  z
    Point3D<ADType> sphere_center; // center of the sphere
    ADType sphere_radius; // radius of the sphere
    auto [cell_x_dim, cell_y_dim, cell_z_dim] = getCellDimensions();
    auto min_cell_dim = std::min(std::min(cell_x_dim, cell_y_dim), cell_z_dim);
    
    ADType rod_radius = min_cell_dim * nCell_rod; // radius of the unchanging rods
    ADType k = k_heaviside/min_cell_dim; // heaviside penalty



      // Redefine global_x_min, global_y_min, and global_z_min
      const double lower_bounding_coord_x = global_x_min - 0.5*min_cell_dim;
      const double lower_bounding_coord_y = global_y_min - 0.5*min_cell_dim;
      const double lower_bounding_coord_z = global_z_min - 0.5*min_cell_dim;

      // get cell dimensions
      auto unit_cell_dim_x = n_cells_into_dual_mesh * min_cell_dim;
      auto unit_cell_dim_y = n_cells_into_dual_mesh * min_cell_dim;
      auto unit_cell_dim_z = n_cells_into_dual_mesh * min_cell_dim;
      
      
      // Loop over cells ... change brGammaSize to number of cells in mesh
      //for (size_t i = 0; i < brGammaSize; i++)
      forAll(mesh.C(),i)
      {

        //Info << "   : Cell : " << i << endl;

        // std::cout << "  : brgamma update!!!\n";
        if (cell_can_change[i]) {
	  //transferVars.brGamma[i] = x[label_index_map[local_cell_id_to_global_var[i]]];

      sphere_radius = unit_cell_dim_x * x[label_index_map[local_cell_id_to_global_var[i]]];

      //if (local_cell_id_to_global_var[i] == 13 ){sphere_radius *= 1.3;}

      /*Info << "   : Cell: " << i 
           << " global label: " << local_cell_id_to_global_var[i]
           << ", sphere_radius: " << sphere_radius[0]
           <<  endl;*/

	  // Get the coordinates of this mesh
	  double xCell = mesh.C()[i].x();
	  double yCell = mesh.C()[i].y();
	  double zCell = mesh.C()[i].z();	  	  
	  
	  Point3D<ADType> point; 
	  point << xCell, yCell, zCell;	  

	  // Set the coodinates of this unit-cell lattice
	  auto label = local_cell_id_to_global_var[i];
	  auto x_int = label / (y_cells * z_cells);
	  auto z_int = (label % (y_cells * z_cells) )/ (y_cells);
	  auto y_int = label % y_cells;

	  // find center of unit-cell
	  double unit_cell_min_x = x_int * unit_cell_dim_x + lower_bounding_coord_x;
	  double unit_cell_min_y = y_int * unit_cell_dim_y + lower_bounding_coord_y;
	  double unit_cell_min_z = z_int * unit_cell_dim_z + lower_bounding_coord_z;	  
	  sphere_center[0] = unit_cell_min_x + unit_cell_dim_x / 2;
	  sphere_center[1] = unit_cell_min_y + unit_cell_dim_y / 2;
	  sphere_center[2] = unit_cell_min_z + unit_cell_dim_z / 2;
	  
	  double sx = unit_cell_min_x + unit_cell_dim_x / 2;
	  double sy = unit_cell_min_y + unit_cell_dim_y / 2;
	  double sz = unit_cell_min_z + unit_cell_dim_z / 2;
	  
	  
	  /*
	  Info << "   : unit_cell_min_x: " << unit_cell_min_x << endl;
	  Info << "   : x_int: " << x_int << endl;
	  Info << "   : unit_cell_dim_x: " << unit_cell_dim_x << endl;	  	  
	  Info << "   : sphere_center[0]: " << unit_cell_min_x + unit_cell_dim_x / 2 << endl;
	  
	  return 0;
	  */
	  
	  //	  sphere_radius = ((gamma_max_value - gamma_min_value) * x[label_index_map[label]] + gamma_min_value) * unit_cell_dim_x / 2.;
	  //sphere_radius = unit_cell_dim_x * 0.3;
	  
	  
	  // set rod points
	  pmx[1] = sphere_center[1]; pmx[2] = sphere_center[2]; pmx[0] = unit_cell_min_x ;
	  ppx[1] = sphere_center[1]; ppx[2] = sphere_center[2]; ppx[0] = unit_cell_min_x + unit_cell_dim_x;

	  pmy[0] = sphere_center[0]; pmy[2] = sphere_center[2]; pmy[1] = unit_cell_min_y;

	  ppy[0] = sphere_center[0]; ppy[2] = sphere_center[2]; ppy[1] = unit_cell_min_y + unit_cell_dim_y;

	  pmz[0] = sphere_center[0]; pmz[1] = sphere_center[1]; pmz[2] = unit_cell_min_z;
	  ppz[0] = sphere_center[0]; ppz[1] = sphere_center[1]; ppz[2] = unit_cell_min_z + unit_cell_dim_z;

	  // evaluate the openfoam cell
	  //gammaIdx[i] = unit_cell_func(point)[0];
	  primitive::Sphere sphere(sphere_radius, sphere_center);
	  primitive::Line rod_x(rod_radius, pmx, ppx);
	  primitive::Line rod_y(rod_radius, pmy, ppy);
	  primitive::Line rod_z(rod_radius, pmz, ppz);
	  auto cross_z = R::Or(R::Or(rod_x, rod_y), rod_z);
	  auto cross_y = R::Or(R::Or(rod_x, rod_z), rod_y);
	  auto cross_x = R::Or(R::Or(rod_y, rod_z), rod_x);
	  auto cross = Average(cross_z, cross_y, cross_x);
	  auto grad_cross = gradient(GenerateFunction(cross), wrt(point), at(point));
	  auto grad_sphere = gradient(GenerateFunction(sphere), wrt(point), at(point));
	  // auto unit_cell = SmoothHeaviside(k, Normalize(cross, grad_cross));
	  // auto unit_cell = SmoothHeaviside(k, Normalize(sphere, grad_sphere));
	  auto unit_cell = SmoothHeaviside(k, R::Or(Normalize(cross, grad_cross), Normalize(sphere, grad_sphere)));
	  auto unit_cell_func = GenerateFunction(unit_cell);

	  ADType f;
	  gammaGeomGrad[i] = -1.0 * unit_cell_dim_x 
                            * gradient(unit_cell_func, wrt(sphere_radius), at(point), f)[0];
	  gammaGeom[i] = 1-f[0];
	  
	  transferVars.brGamma[i] = gammaGeom[i];
	  //transferVars.brGammaStart[i] = gammaGeom[i];	  don't think this necessary
    
        }
    }
      

   
    gammaGeomGrad.correctBoundaryConditions();
    gammaGeom.correctBoundaryConditions();
      
    newDesignParameters();
    evalSensitivities(); // remove for final code
    evalObjective(); //remove for final code
    saveResults();
    


    };



    opt.UpdatedVariableCallback();
    




    // Call once to initialize problem
    if (no_design_vars)
    {
      std::cout << "  : OP pre update variable callback; No design vars on rank: " << my_rank << "\n";
    }
 
    if (no_design_vars)
    {
      std::cout << "  : Op variable callback works; No design vars on rank: " << my_rank << "\n";
    }
    // method we'll call to go
    opt.go.onPreprocess([&]() {
        // set objective
        nlopt_options.Double["constraint_tol"] = 1.e-8;
        opt.setObjective(obj);
        nlopt_options.Double["constraint_tol"] = 1.e-8;
        opt.addConstraint(cons1);
      });
    if (no_design_vars)
    {
      std::cout << "  : op go called; No design vars on rank: " << my_rank << "\n";
    }
    // Run optimizer
    try {
      opt.Go();
      std::cout << "found minimum = " << std::setprecision(10) << opt.Solution() / scale
                << std::endl;
      // make sure to save openfoam results of the solution
      opt.UpdatedVariableCallback();
      saveResults();
    } catch (std::exception& e) {
      std::cout << "nlopt failed: " << e.what() << std::endl;
    }

    delete [] transferVars.brGammaStart;
    delete [] transferVars.brGamma;
    delete [] transferVars.brSensitivities;
    delete [] transferVars.brUpperBound;
    delete [] transferVars.brLowerBound;
    delete [] transferVars.VCsens;
    return 0;
    
}

void saveResults()
{
    Foam::Time              & myRunTime         = *transferVars.runTime;
    volScalarField          & gamma             = *transferVars.gamma;
    volScalarField          & gammaRaw             = *transferVars.gammaRaw;
    double                  & brOptimizationObjective = *transferVars.brOptimizationObjective;       
    double                  & VC                = *transferVars.VC;       
    volScalarField          & gammaId             = *transferVars.gammaId;
    volScalarField          & cellChangeMask             = *transferVars.cellChangeMask;
    volScalarField          & gammaIdx             = *transferVars.gammaIdx;
	  volScalarField          & gammaIdz             = *transferVars.gammaIdz;
    volScalarField          & rank             = *transferVars.rank;
 	  volScalarField          & rank_label             = *transferVars.rank_label;

	  ++myRunTime;
    Info << nl << "\n\n\n  : Time New = " << myRunTime.timeName() << " objective: " << brOptimizationObjective << " constraint: " << VC <<endl;
    //gamma.write();  // only write gamma
    //gammaRaw.write();
    // The following volScalarFields are useful for debugging
    //gammaId.write();
    //gammaIdx.write();
    //gammaIdz.write();
    //cellChangeMask.write();
    //rank.write();
    //rank_label.write();
    // dump();
    // gammaRaw.write();
    myRunTime.write();
}

void saveResults(std::function<void()> & dump)
{
    Foam::Time              & myRunTime         = *transferVars.runTime;
    volScalarField          & gamma             = *transferVars.gamma;
    volScalarField          & gammaRaw             = *transferVars.gammaRaw;
    double                  & brOptimizationObjective = *transferVars.brOptimizationObjective;       
    double                  & VC                = *transferVars.VC;       
    // Info << nl << "  : Time = " << myRunTime.timeName() << nl << endl;
    ++myRunTime;
    Info << nl << "  : Time New = " << myRunTime.timeName() << " objective: " << brOptimizationObjective << " constraint: " << VC <<endl;
    //gamma.write();  // only write gamma
    dump();
    // gammaRaw.write();
    myRunTime.write();
}

double evalConstraint()
{

  volScalarField          & gamma             = *transferVars.gamma;
  const volScalarField::Internal & meshVolumes       = *transferVars.meshVolumes;
  double                  & VC                = *transferVars.VC;       

  VC = 1.0 - gSum(gamma*1.0*meshVolumes) / gSum(meshVolumes);
  std::cout << "  : Constraint: " << VC << " \n";

  return VC;
}

void evalConstraintSensitivity()
{
  const volScalarField::Internal & meshVolumes       = *transferVars.meshVolumes;
  double                  & VC                = *transferVars.VC;
  double                  & VT                = *transferVars.VT;
  const cellList          & cells             = *transferVars.cells;
  volScalarField          & gammaGeomGrad             = *transferVars.gammaGeomGrad;
  dimensionedScalar       & min_gamma         = *transferVars.min_gamma;                                 
  dimensionedScalar       & max_gamma         = *transferVars.max_gamma;

  double vc_vt = -1.0*(max_gamma.value()-min_gamma.value())/gSum(meshVolumes);

  forAll(cells, cellI)
  {
    transferVars.VCsens[cellI] = vc_vt*meshVolumes[cellI]*gammaGeomGrad[cellI]; 
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

        
    
    
    brOptimizationObjective = gSum((topI2 & topI2.patch().Sf())
                        *(-1.*topPhi2 - stdPotential.value()) )
                        + gSum(-1000.*inletPressure*inPhi);
    
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
    double                  & brGammaMin        = *transferVars.brGammaMin;
    double                  & brGammaMax        = *transferVars.brGammaMax;
    auto                    & cell_can_change   = *transferVars.cell_can_change;
    auto                    & is_top_cell       = *transferVars.is_top_cell;
    dimensionedScalar       & min_gamma         = *transferVars.min_gamma;                                 
    dimensionedScalar       & max_gamma         = *transferVars.max_gamma;

  forAll(cells, cellI)
  {
    if (cell_can_change[cellI]) {
      transferVars.brGammaStart[cellI] = starting_value; 
      transferVars.brGamma[cellI] = starting_value;
    }
    else {
      if (is_top_cell[cellI]) {
        transferVars.brGammaStart[cellI] = min_gamma.value();
        transferVars.brGamma[cellI] = min_gamma.value();
      }
      else {
        transferVars.brGammaStart[cellI] = max_gamma.value();
        transferVars.brGamma[cellI] = max_gamma.value();
      }
    }
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
    volScalarField          & gammaAux             = *transferVars.gammaAux;
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
    double                  & qArea             = *transferVars.qArea;
    double                  & qPerm             = *transferVars.qPerm;    
    dimensionedScalar       & a0               = *transferVars.a0;
    dimensionedScalar       & porosity0               = *transferVars.porosity0;
    dimensionedScalar       & permeability0               = *transferVars.permeability0;    
    dimensionedScalar       & stdPotential      = *transferVars.stdPotential;
    volScalarField          & gammaGeomGrad             = *transferVars.gammaGeomGrad;

    dimensionedScalar correctDim2("correctDim2",
                                  dimensionSet(0, -2, 3, 0, 0, 0, 0),1.0);        
        
             
    deltaPhi = Phi1 -Phi2 - stdPotential;
            
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


        sensPhi =   -1.0 * a0 * qArea 
                         * ( qArea + 1.0) * ( max_gamma - min_gamma)
                         / pow(   
                                qArea * ( max_gamma - min_gamma) 
                                + gamma - min_gamma
                              ,2) 
                    * exchangeCurrent * (Psi1 - Psi2)                    
                    - 1.5 * (1.-porosity0) 
                          * sqrt( (1. - porosity0) * (1. -gamma) )* kappa1
                          * (fvc::grad(Phi1)&fvc::grad(Psi1))
                    + 1.5 * (1.-porosity0) 
                          * sqrt( gamma + porosity0 * (1. - gamma) ) * kappa2                    
                          * (fvc::grad(Phi2)&fvc::grad(Psi2));

                    
        
        sensConc =  -1.0 * a0 * qArea 
                         * ( qArea + 1.0) * ( max_gamma - min_gamma)
                         / pow(   
                                qArea * ( max_gamma - min_gamma) 
                                + gamma - min_gamma
                              ,2) 
                    * exchangeCurrent * (CRa - COa)/nelec/faradayConstant
                    + 1.5 * (1.-porosity0) 
                          * sqrt( gamma + porosity0 * (1. - gamma) ) * DCO                    
                          * (fvc::grad(CO)&fvc::grad(COa))*CRef
                    + 1.5 * (1.-porosity0) 
                          * sqrt( gamma + porosity0 * (1. - gamma) ) * DCR                    
                          * (fvc::grad(CR)&fvc::grad(CRa))*CRef;        


        sensVel = -0.5 * nu/permeability0 * (U&Ua)
            * betaTanh / (max_gamma-min_gamma) / Foam::tanh(betaTanh*0.5) 
            * ( 1.0 - Foam::pow( 
                        Foam::tanh(
                            betaTanh*(   (gamma-min_gamma)
                                       / (max_gamma-min_gamma) - 0.5)
                                  )
                        , 2) 
              );

        
        sens = sensPhi + sensConc + sensVel*correctDim2;


        sensRaw = sens * (max_gamma-min_gamma);
        
    sensRaw.correctBoundaryConditions();

    

    
    forAll(cells, cellI)
    {
        sensRaw[cellI]=meshVol[cellI]*sensRaw[cellI]*gammaGeomGrad[cellI];
        if (centroids[cellI][1] < 0.0)
        {
            sensRaw[cellI]=0.0;
        }
    }
    


    forAll(cells, cellI)
    {
        transferVars.brSensitivities[cellI] = sensRaw[cellI];
    }
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
    dimensionedScalar          & bR                = *transferVars.bR;
    dimensionedScalar          & bO                = *transferVars.bO;
    dimensionedScalar          & gR                = *transferVars.gR;
    dimensionedScalar          & gO                = *transferVars.gO;   
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
    auto                    & cell_can_change      = *transferVars.cell_can_change;
    double                  & qArea             = *transferVars.qArea;
    double                  & qPerm             = *transferVars.qPerm;    
    dimensionedScalar       & a0               = *transferVars.a0;
    dimensionedScalar       & porosity0               = *transferVars.porosity0;
    dimensionedScalar       & permeability0               = *transferVars.permeability0;    
    volScalarField          & exchangeCurrent   = *transferVars.exchangeCurrent;

    /// First, update design variables with output of optimizer
    // gammaRaw is expected to be in [0,1]
    
    forAll(cells, cellI)
    {
        gammaRaw[cellI] = transferVars.brGamma[cellI];
    }

    gammaRaw.correctBoundaryConditions();
 
    
    
    // Regularize the updated design variable (gammaRaw)     
    
    gammaAux = gammaRaw;        
    gammaAux = min(max(gammaAux, 0.0),1.0);
    gamma = (max_gamma - min_gamma) * gammaAux + min_gamma;                         
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

    Info << "Max_gamma: " << max_gamma << endl;
    Info << "Min_gamma: " << min_gamma << endl;    
    Info << "qPerm: " << qPerm << endl;        
    Info << "qArea: " << qArea << endl;            
    
    
    
    porosity = gamma+porosity0*(1.0-gamma);
         

    invPermeability = 0.5/permeability0  *
            (  1.0 
              - Foam::tanh(betaTanh * (
                                        (gamma - min_gamma)/(max_gamma - min_gamma)
                                        - 0.5
                                      )
                          )
            / ( Foam::tanh(0.5*betaTanh) ));   
                           
    
     areaPerVol = a0* qArea*(max_gamma - gamma)
                       /
                        (qArea*(max_gamma - min_gamma)+gamma-min_gamma);

    
    forAll(cells, cellI)
    {
        
            
    
            if (centroids[cellI][1] < 0.0)
            {
                porosity[cellI] = 1.0 - 1.0e-4;
                areaPerVol[cellI] = 0.0;
                invPermeability[cellI] = 0.0;
            }
            // invPermeability should be read in
            if (centroids[cellI][1] < -0.0012)
            {
                 invPermeability[cellI] = 0.0;
            }

            
        

    }
    
    porosity.correctBoundaryConditions();
    
    

        
        permeability = 1.0 / (invPermeability + 1.0e-20/permeabilityRef);        
        effDCO=sqrt(Foam::pow(porosity,3))*DCO;
        effDCR=sqrt(Foam::pow(porosity,3))*DCR;
        effKappa1=sqrt(Foam::pow(1-porosity,3))*kappa1;
        effKappa2=sqrt(Foam::pow(porosity,3))*kappa2;   
        
      
                
        int iterVel = 0;
        double flowConsErr = 1.0;        
                                
        // Pressure-velocity SIMPLE corrector

       
        
        double flowRes = 1.1;
	      while ( (iterVel < interation_limit && flowRes > solver_tol) || iterVel < 5)
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
         
            Info << "gR: " << gR << endl;
            Info << "gO: " << gO << endl;            
            Info << "bR: " << bR << endl;
            Info << "bO: " << bO<< endl;      
            Info << "velRef: " << velRef << endl;     
 

         
            kmR = gR * Foam::pow(mag(U/velRef)+1.e-20,bR)+1.0e-20*velRef;
            kmO = gO * Foam::pow(mag(U/velRef)+1.e-20,bO)+1.0e-20*velRef;           
                       
            kmR.correctBoundaryConditions();
            kmO.correctBoundaryConditions();


            Info << kmO << endl;

            int iterPot = 0;
            double currConsErr = 1.0; 
            double concConsErr = 1.0;   
            
            double maxScalarsRes = 1.0;
            while ( (iterPot < interation_limit && maxScalarsRes > solver_tol) || iterPot < 5 )
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
            
            while ( (iterAdjnts < interation_limit && maxAdjRes > solver_tol) || iterAdjnts < 5)
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
        
        while ( (iterVelAdjnts < interation_limit && maxVelAdjRes > solver_tol) || iterVelAdjnts < 5)     
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
