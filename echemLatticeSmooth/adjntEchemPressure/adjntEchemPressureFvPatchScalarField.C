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

#include "adjntEchemPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "turbulentTransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjntEchemPressureFvPatchScalarField::
adjntEchemPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


Foam::adjntEchemPressureFvPatchScalarField::
adjntEchemPressureFvPatchScalarField
(
    const adjntEchemPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::adjntEchemPressureFvPatchScalarField::
adjntEchemPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF)
{
    fvPatchField<scalar>::operator=
    (
        scalarField("value", dict, p.size())
    );
}


Foam::adjntEchemPressureFvPatchScalarField::
adjntEchemPressureFvPatchScalarField
(
    const adjntEchemPressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::adjntEchemPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>("phi");

    const fvsPatchField<scalar>& phiap =
        patch().lookupPatchField<surfaceScalarField, scalar>("phia");

    const fvPatchField<vector>& Up =
        patch().lookupPatchField<volVectorField, vector>("U");

    const fvPatchField<vector>& Uap =
        patch().lookupPatchField<volVectorField, vector>("Ua");


    // Normal of primaly velocity. Take momentum flux, phip = u*Area, 
    // at the patch surface and divide by area

    
    scalarField Up_n = phip / patch().magSf(); //primal
    
    scalarField Uap_n = phiap / patch().magSf(); //primal 

    // Uaneigh_n is the node neighboring the external path. nf = face normal.
    scalarField Uaneigh_n = (Uap.patchInternalField() & patch().nf());


// Need to lookup the viscosity... infield form.
    // Extract the dictionary from the database
    const dictionary& transportProperties = db().lookupObject<IOdictionary>
    (
        "transportProperties"
    );
     // Extracting scalar value
    const dimensionedScalar nu(transportProperties.lookup("nu"));
    scalarField nueff(Uaneigh_n.size(),nu.value());    
    const scalarField& deltainv = patch().deltaCoeffs();
    




    operator==(
        //(Up & Uap)  //important to keep dot product in paranthesis 
          Up_n * Uap_n 
        + nueff * deltainv * (Uap_n - Uaneigh_n)
        //+ C1a*C1
        );

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::adjntEchemPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        adjntEchemPressureFvPatchScalarField
    );
}

// ************************************************************************* //
