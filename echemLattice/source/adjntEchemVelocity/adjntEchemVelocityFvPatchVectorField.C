/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "adjntEchemVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvPatchFieldMapper.H"
#include "turbulentTransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjntEchemVelocityFvPatchVectorField::
adjntEchemVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


Foam::adjntEchemVelocityFvPatchVectorField::
adjntEchemVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF)
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::adjntEchemVelocityFvPatchVectorField::
adjntEchemVelocityFvPatchVectorField
(
    const adjntEchemVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::adjntEchemVelocityFvPatchVectorField::
adjntEchemVelocityFvPatchVectorField
(
    const adjntEchemVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void Foam::adjntEchemVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<scalar>& phiap =
        patch().lookupPatchField<surfaceScalarField, scalar>("phia");

    const fvPatchField<vector>& Up =
        patch().lookupPatchField<volVectorField, vector>("U");
        
    const fvPatchField<vector>& Uap =
        patch().lookupPatchField<volVectorField, vector>("Ua");
        
    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>("phi");
        
    
    const scalarField& deltainv = patch().deltaCoeffs();
    
    // Normal of primaly velocity. Take momentum flux, phip = u*Area, 
    // at the patch surface and divide by area

    scalarField Up_ns = phip / patch().magSf(); //primal
        
    vectorField Uaneigh = Uap.patchInternalField();

    // Need to lookup the viscosity... infield form.
    // Extract the dictionary from the database
    const dictionary& transportProperties = db().lookupObject<IOdictionary>
    (
        "transportProperties"
    );
     // Extracting scalar value
    const dimensionedScalar nu(transportProperties.lookup("nu"));
    scalarField nueff(Uaneigh.size(),nu.value());   
    
    vectorField Uaneigh_n = (Uaneigh & patch().nf()) * patch().nf();
    
    vectorField Uaneigh_t = Uaneigh - Uaneigh_n;
    
    vectorField Uap_t = (nueff*deltainv*Uaneigh_t)
                        / (Up_ns + nueff*deltainv);

    vectorField Uap_n = (phiap * patch().Sf())
                       / (patch().magSf()*patch().magSf());
    

    vectorField::operator==(
    
        Uap_t + Uap_n
    
    );
    //vectorField::operator=(Uan + UtHat);

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::adjntEchemVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        adjntEchemVelocityFvPatchVectorField
    );
}


// ************************************************************************* //
