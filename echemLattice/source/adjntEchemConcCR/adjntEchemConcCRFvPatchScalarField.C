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

#include "adjntEchemConcCRFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::adjntEchemConcCRFvPatchScalarField::t() const
{
    return db().time().timeOutputValue();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjntEchemConcCRFvPatchScalarField::
adjntEchemConcCRFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF)
{
    //refGrad() = 0.0;
    //refValue() = 4.2;
    //valueFraction() = 0.99999;
}


Foam::adjntEchemConcCRFvPatchScalarField::
adjntEchemConcCRFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF)
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


Foam::adjntEchemConcCRFvPatchScalarField::
adjntEchemConcCRFvPatchScalarField
(
    const adjntEchemConcCRFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::adjntEchemConcCRFvPatchScalarField::
adjntEchemConcCRFvPatchScalarField
(
    const adjntEchemConcCRFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf)
{}


Foam::adjntEchemConcCRFvPatchScalarField::
adjntEchemConcCRFvPatchScalarField
(
    const adjntEchemConcCRFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::adjntEchemConcCRFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
}


void Foam::adjntEchemConcCRFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const adjntEchemConcCRFvPatchScalarField& tiptf =
        refCast<const adjntEchemConcCRFvPatchScalarField>(ptf);
}


void Foam::adjntEchemConcCRFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField& CR =
        patch().lookupPatchField<volScalarField, scalar>("CR");
    const fvsPatchScalarField& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>("phi");
    const fvPatchField& porosity =
        patch().lookupPatchField<volScalarField, scalar>("porosity");
    
    scalarField vp_n = phip / this->patch().magSf(); //primal
    

    // Need to lookup the diffusivity... in field form.
    // Extract the dictionary from the database
    const dictionary& transportProperties = db().lookupObject<IOdictionary>
    (
        "transportProperties"
    );
    const dimensionedScalar DCR(transportProperties.lookup("DCR"));
    scalarField DCR_field(vp_n.size(),DCR.value());

    DCR_field = DCR_field*sqrt(pow(porosity,3));

    scalarField W = vp_n/DCR_field/this->patch().deltaCoeffs();



    refGrad() = 0.0;   
    refValue() = 0.0;
    valueFraction() = W / (1.0 + W) ;

    
    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::adjntEchemConcCRFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        adjntEchemConcCRFvPatchScalarField
    );
}

// ************************************************************************* //
