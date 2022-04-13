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

#include "adjntCOObjCOFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::adjntCOObjCOFvPatchScalarField::t() const
{
    return db().time().timeOutputValue();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjntCOObjCOFvPatchScalarField::
adjntCOObjCOFvPatchScalarField
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


Foam::adjntCOObjCOFvPatchScalarField::
adjntCOObjCOFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    TT_(readScalar(dict.lookup("targetConc")))
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


Foam::adjntCOObjCOFvPatchScalarField::
adjntCOObjCOFvPatchScalarField
(
    const adjntCOObjCOFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    TT_(ptf.TT_)
{}


Foam::adjntCOObjCOFvPatchScalarField::
adjntCOObjCOFvPatchScalarField
(
    const adjntCOObjCOFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf)
{}


Foam::adjntCOObjCOFvPatchScalarField::
adjntCOObjCOFvPatchScalarField
(
    const adjntCOObjCOFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    TT_(ptf.TT_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::adjntCOObjCOFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
}


void Foam::adjntCOObjCOFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const adjntCOObjCOFvPatchScalarField& tiptf =
        refCast<const adjntCOObjCOFvPatchScalarField>(ptf);
}


void Foam::adjntCOObjCOFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    
    const fvPatchField& CO =
        patch().lookupPatchField<volScalarField, scalar>("CO");
    /*const fvsPatchScalarField& phip =
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
    const dimensionedScalar DCO(transportProperties.lookup("DCO"));
    scalarField DCO_field(vp_n.size(),DCO.value());

    DCO_field = DCO_field*sqrt(pow(porosity,3));

    scalarField W = vp_n/DCO_field/this->patch().deltaCoeffs();*/


    const fvPatchField& DCO_field =
        patch().lookupPatchField<volScalarField, scalar>("effDCO");

    refGrad() = (TT_ - CO) / DCO_field;   
    refValue() = 0.0;
    valueFraction() = 0.0 ;

    
    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::adjntCOObjCOFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("targetConc") << TT_ << token::END_STATEMENT << nl;    
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        adjntCOObjCOFvPatchScalarField
    );
}

// ************************************************************************* //
