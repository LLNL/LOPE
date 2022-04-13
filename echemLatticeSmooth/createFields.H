Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

Info<< "Reading field Ua\n" << endl;
volVectorField Ua
(
    IOobject
    (
        "Ua",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading/calculating face flux field phia\n" << endl;

surfaceScalarField phia
(
    IOobject
    (
        "phia",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(Ua) & mesh.Sf()
);

Info<< "Reading field pa\n" << endl;
volScalarField pa
(
    IOobject
    (
        "pa",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, simple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

label paRefCell = 0;
scalar paRefValue = 0.0;
setRefCell
(
    pa,
    simple.dict(),
    paRefCell,
    paRefValue
);
mesh.setFluxRequired(pa.name());

Info<< "Reading field CO\n" << endl;
volScalarField CO
(
    IOobject
    (
        "CO",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field CR\n" << endl;
volScalarField CR
(
    IOobject
    (
        "CR",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field COa\n" << endl;
volScalarField COa
(
    IOobject
    (
        "COa",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field CRa\n" << endl;
volScalarField CRa
(
    IOobject
    (
        "CRa",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
    
Info<< "Reading field Phi1\n" << endl;
volScalarField Phi1
(
    IOobject
    (
        "Phi1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Phi2\n" << endl;
volScalarField Phi2
(
    IOobject
    (
        "Phi2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Psi1\n" << endl;
volScalarField Psi1
(
    IOobject
    (
        "Psi1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
    
Info<< "Reading field Psi2\n" << endl;
volScalarField Psi2
(
    IOobject
    (
        "Psi2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

    
dimensionedScalar nu
(
    transportProperties.lookup("nu")
);

dimensionedScalar DCO
(
    transportProperties.lookup("DCO")
);
    
dimensionedScalar DCR
(
    transportProperties.lookup("DCR")
);
       
dimensionedScalar kappa1
(
    transportProperties.lookup("kappa1")
);

dimensionedScalar kappa2
(
    transportProperties.lookup("kappa2")
);

dimensionedScalar i0
(
    transportProperties.lookup("i0")
);

dimensionedScalar Temperature
(
    transportProperties.lookup("Temperature")
);
    
dimensionedScalar stdPotential
(
    transportProperties.lookup("stdPotential")
);
    
dimensionedScalar alphaA
(
    transportProperties.lookup("alphaA")
);
    
dimensionedScalar alphaC
(
    transportProperties.lookup("alphaC")
);

dimensionedScalar nelec
(
    transportProperties.lookup("nelec")
);

dimensionedScalar CRef
(
    transportProperties.lookup("CRef")
);
  
dimensionedScalar a0
(
    transportProperties.lookup("a0")
);

dimensionedScalar porosity0
(
    transportProperties.lookup("porosity0")
);


dimensionedScalar rFiber
(
    transportProperties.lookup("rFiber")
);

/*
dimensionedScalar gO
(
    transportProperties.lookup("gO")
);

dimensionedScalar gR
(
    transportProperties.lookup("gR")
);

dimensionedScalar bO
(
    transportProperties.lookup("bO")
);

dimensionedScalar bR
(
    transportProperties.lookup("bR")
);
*/

dimensionedScalar unitCellLen
(
    transportProperties.lookup("unitCellLen")
);

dimensionedScalar printResLen
(
    transportProperties.lookup("printResLen")
);

dimensionedScalar stepLength
(
    transportProperties.lookup("stepLength")
);

dimensionedScalar stepLengthInit
(
    transportProperties.lookup("stepLengthInit")
);  

dimensionedScalar stepLengthMin
(
    transportProperties.lookup("stepLengthMin")
);  

dimensionedScalar stepLengthIter
(
    transportProperties.lookup("stepLengthIter")
);

dimensionedScalar gradPenalty
(
    transportProperties.lookup("gradPenalty")
);

dimensionedScalar velocityPenalty
(
    transportProperties.lookup("velocityPenalty")
);

dimensionedScalar rFilter
(
    transportProperties.lookup("rFilter")
);

dimensionedScalar volPenalty
(
    transportProperties.lookup("volPenalty")
);

dimensionedScalar qPermInterp0_
(
    transportProperties.lookup("qPermInterp0")
);

dimensionedScalar qPermInterpFinal_
(
    transportProperties.lookup("qPermInterpFinal")
);

dimensionedScalar NqPermInterp_
(
    transportProperties.lookup("NqPermInterp")
);

dimensionedScalar qArea0_
(
    transportProperties.lookup("qArea0")
);

dimensionedScalar qAreaFinal_
(
    transportProperties.lookup("qAreaFinal")
);

dimensionedScalar NqArea_
(
    transportProperties.lookup("NqArea")
);

dimensionedScalar betaTanh0_
(
    transportProperties.lookup("betaTanh0")
);

dimensionedScalar betaTanhMax_
(
    transportProperties.lookup("betaTanhMax")
);

dimensionedScalar NbetaTanh_
(
    transportProperties.lookup("NbetaTanh")
);

dimensionedScalar etaTanh
(
    transportProperties.lookup("etaTanh")
);

dimensionedScalar gammaTrnPermInterp0_
(
    transportProperties.lookup("gammaTrnPermInterp0")
);

dimensionedScalar gammaTrnPermInterpFinal_
(
    transportProperties.lookup("gammaTrnPermInterpFinal")
);

dimensionedScalar deltaNSteps_
(
    transportProperties.lookup("deltaNSteps")
);

dimensionedScalar permeability0
(
    transportProperties.lookup("permeability0")
);
      
dimensionedScalar gasConstant
("gasConstant", dimensionSet(1, 2, -2, -1, 0, 0, 0), 
    8.3144598
);
    
dimensionedScalar faradayConstant
("faradayConstant", dimensionSet(0, 0, 1, 0, 0, 1, 0), 
    96485.3329
);
    
dimensionedScalar FRT 
("FRT", faradayConstant/Temperature/gasConstant);

volScalarField gamma
(
    IOobject
    (
       "gamma",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
  
volScalarField porosity
(
    IOobject
    (
       "porosity",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
   
volScalarField areaPerVol
(
    IOobject
    (
        "areaPerVol",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    a0*gamma
);

volScalarField invPermeability
(
    IOobject
    (
        "invPermeability",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    1.0/permeability0 * gamma
);

dimensionedScalar permeabilityRef
(
    "permeabilityRef", 
    dimensionSet(0, 2, 0, 0, 0, 0, 0),1.0
);

volScalarField permeability
(
    IOobject
    (
        "permeability",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    1.0 / (invPermeability + 1.0e-20/permeabilityRef)
);

const dimensionedScalar velRef
(
    "velRef", 
    dimensionSet(0, 1, -1, 0, 0, 0, 0),1.0
);

const dimensionedScalar lenRef
(
    "lenRef", 
    dimensionSet(0, 1, 0, 0, 0, 0, 0),1.0
);

const dimensionedScalar c1_km("c1_km", dimless, 0.468);
const dimensionedScalar c2_km("c2_km", dimless, 0.391);
const dimensionedScalar c3_km("c3_km", dimless, 0.4672);
const dimensionedScalar c4_km("c4_km", dimless, -0.1247);

volScalarField bR
(
    IOobject
    (
        "bR",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    c4_km*porosity
);

volScalarField bO
(
    IOobject
    (
        "bO",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    c4_km*porosity
);


volScalarField kmR
(
    IOobject
    (
        "kmR",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);    
      
volScalarField kmO
(
    IOobject
    (
        "kmO",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);          

volScalarField gR
(
    IOobject
    (
        "gR",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    kmO
);

volScalarField gO
(
    IOobject
    (
        "gO",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    kmR
);
        
volScalarField deltaPhi
(
    IOobject
    (
        "deltaPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    Phi1 -Phi2 - stdPotential
);    
    
Info<< "Creating field I2\n" << endl;
volVectorField I2
(
    IOobject
    (
        "I2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    -1.0*kappa2*fvc::grad(Phi1)
);

Info<< "Creating field I1\n" << endl;
volVectorField I1
(
    IOobject
    (
        "I1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    -1.0*kappa1*fvc::grad(Phi1)
);

volScalarField exchangeCurrent
(
    IOobject
    (
        "exchangeCurrent",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    i0 * (CR - CO + FRT * deltaPhi * (CR*alphaA + CO*alphaC))
);
  
volScalarField effDCO
(
    IOobject
    (
        "effDCO",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    sqrt(pow(porosity,3))*DCO
);  

volScalarField effDCR
(
    IOobject
    (
        "effDCR",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    sqrt(pow(porosity,3))*DCR
);

volScalarField effKappa1
(
    IOobject
    (
        "effKappa1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    sqrt(pow(1-porosity,3))*kappa1
);

volScalarField effKappa2
(
    IOobject
    (
        "effKappa2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    sqrt(pow(porosity,3))*kappa2
); 
    
volScalarField sensPhi
(
    IOobject
    (
        "sensPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    a0 * exchangeCurrent * (Psi1 - Psi2)
);

volScalarField sensConc
(
    IOobject
    (
        "sensConc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    a0*exchangeCurrent*(CRa - COa)/nelec/faradayConstant
);

volScalarField sensVel
(
    IOobject
    (
        "sensVel",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (U&Ua) * 24. * nu / sqr(rFiber) * sqrt(pow(1.-porosity0,3)) * sqrt(gamma)
           * (1. + 168. * pow(gamma*(1.-porosity0),3))               
);

volScalarField sens
(
    IOobject
    (
        "sens",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    sensPhi + sensConc
);

volScalarField sensRaw
(
    IOobject
    (
        "sensRaw",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField gammaAux
(
    IOobject
    (
        "gammaAux",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField gammaRaw
(
    IOobject
    (
        "gammaRaw",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    gamma
);

volScalarField gammaOld
(
    IOobject
    (
        "gammaOld",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    gamma
);

volScalarField gammaRawOld
(
    IOobject
    (
        "gammaRawOld",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    gamma
);

volScalarField deltaGamma
(
    IOobject
    (
        "deltaGamma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    gamma
);

volScalarField rFiberCell
(
    IOobject
    (
        "rFiberCell",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    gamma*rFiber
);

volScalarField dPorosity
(
    IOobject
    (
        "dPorosity",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    porosity/rFiber
);

volScalarField dAreaPerVol
(
    IOobject
    (
        "dAreaPerVol",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    areaPerVol/rFiber
);

volScalarField dInvPermeability
(
    IOobject
    (
        "dInvPermeability",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    invPermeability/rFiber
);


labelList inletCells = mesh.boundary()["inlet"].faceCells();
labelList outletCells = mesh.boundary()["outlet"].faceCells();
labelList topCells = mesh.boundary()["top"].faceCells();
labelList btmCells = mesh.boundary()["btm"].faceCells();



