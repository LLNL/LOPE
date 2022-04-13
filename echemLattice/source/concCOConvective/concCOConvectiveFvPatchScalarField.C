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

#include "concCOConvectiveFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::concCOConvectiveFvPatchScalarField::t() const
{
    return db().time().timeOutputValue();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::concCOConvectiveFvPatchScalarField::
concCOConvectiveFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    TT_(0.0),
    PP_(0.0)
{
    //refGrad() = 0.0;
    //refValue() = 4.2;
    //valueFraction() = 0.99999;
}


Foam::concCOConvectiveFvPatchScalarField::
concCOConvectiveFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    TT_(readScalar(dict.lookup("farFieldConc"))),
    PP_(readScalar(dict.lookup("massTransferCoef")))    
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


Foam::concCOConvectiveFvPatchScalarField::
concCOConvectiveFvPatchScalarField
(
    const concCOConvectiveFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    TT_(ptf.TT_),
    PP_(ptf.PP_)    
{}


Foam::concCOConvectiveFvPatchScalarField::
concCOConvectiveFvPatchScalarField
(
    const concCOConvectiveFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf)
{}


Foam::concCOConvectiveFvPatchScalarField::
concCOConvectiveFvPatchScalarField
(
    const concCOConvectiveFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    TT_(ptf.TT_),
    PP_(ptf.PP_)    
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::concCOConvectiveFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
}


void Foam::concCOConvectiveFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const concCOConvectiveFvPatchScalarField& tiptf =
        refCast<const concCOConvectiveFvPatchScalarField>(ptf);
}


void Foam::concCOConvectiveFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    /*
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
    */
    
    const fvPatchField& DCO_field =
        patch().lookupPatchField<volScalarField, scalar>("effDCO");


    refGrad() = 0.0 ;   
    refValue() = TT_;
    valueFraction() = 1./(1.+DCO_field*this->patch().deltaCoeffs()/PP_);

    
    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::concCOConvectiveFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("farFieldConc") << TT_ << token::END_STATEMENT << nl;
    os.writeKeyword("massTransferCoef") << PP_ << token::END_STATEMENT << nl;    
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        concCOConvectiveFvPatchScalarField
    );
}

// ************************************************************************* //
