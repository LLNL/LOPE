/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

#include "fixedCurrentFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::fixedCurrentFvPatchScalarField::t() const
{
    return db().time().timeOutputValue();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedCurrentFvPatchScalarField::
fixedCurrentFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    TT_(0.0)
{
    //refGrad() = 0.0;
    //refValue() = 4.2;
    //valueFraction() = 0.99999;
}


Foam::fixedCurrentFvPatchScalarField::
fixedCurrentFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    TT_(readScalar(dict.lookup("currentFlux")))
{
    //refGrad() = Zero;
    //valueFraction() = 0.0;

    //refValue() = scalarField("value", dict, p.size());
        
    //fvPatchScalarField::operator=(patchInternalField());

    /*
    fvPatchScalarField::operator=
    (
        //patchInternalField()
        valueFraction()*refValue()
      +
        (1.0 - valueFraction())*
        (
            this->patchInternalField()
          + refGrad()/this->patch().deltaCoeffs()
        )
    );
    */

    //evaluate();
    
    //evaluate();

         

    
    //Initialise with the value entry if evaluation is not possible
    fvPatchScalarField::operator=
    (
        scalarField("value", dict, p.size())
    );
    this->refValue() = *this;
    
}


Foam::fixedCurrentFvPatchScalarField::
fixedCurrentFvPatchScalarField
(
    const fixedCurrentFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    TT_(ptf.TT_)
{}


Foam::fixedCurrentFvPatchScalarField::
fixedCurrentFvPatchScalarField
(
    const fixedCurrentFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf)
{}


Foam::fixedCurrentFvPatchScalarField::
fixedCurrentFvPatchScalarField
(
    const fixedCurrentFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    TT_(ptf.TT_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedCurrentFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
}


void Foam::fixedCurrentFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const fixedCurrentFvPatchScalarField& tiptf =
        refCast<const fixedCurrentFvPatchScalarField>(ptf);
}


void Foam::fixedCurrentFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField& Phi2 =
        patch().lookupPatchField<volScalarField, scalar>("Phi2");
   
    const fvPatchField& porosity =
        patch().lookupPatchField<volScalarField, scalar>("porosity");
   
    // Need to lookup the diffusivity... in field form.
    // Extract the dictionary from the database
    const dictionary& transportProperties = db().lookupObject<IOdictionary>
    (
        "transportProperties"
    );
    const dimensionedScalar kappa2(transportProperties.lookup("kappa2"));

    scalarField kappa2_field(Phi2.size(),kappa2.value());
    
    kappa2_field = kappa2_field*sqrt(pow(porosity,3));
    //scalarField W = vp_n/DC1_field/this->patch().deltaCoeffs();



    refGrad() = TT_ / kappa2_field ;   
    refValue() = 0.0;
    valueFraction() = 0.0;

    
    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::fixedCurrentFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("currentFlux") << TT_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedCurrentFvPatchScalarField
    );
}

// ************************************************************************* //
