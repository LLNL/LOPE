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

Class
    Foam::adjntCOObjCOFvPatchScalarField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides an external convection BC for Laplaces 
    equation assuming NO SOURCE. The reason we require NO SOURCE is that this
    implementation is still incomplete. The conductivity should be specified in
    the transportProperties dictionary and not in the boundary condition. 
    
    However, if we assume NO SOURCE, then the conductivity appears only in the 
    BC. That is what is implemented here.
    
    Calculated as:

        \f[
            -k dT/dn = h * (T - Tinf)
        \f]

    where
    \vartable
        k      | conductivity [units]
        htc    | heat transfer coefficient [units]
        Tinf   | exterior temperature [units]
    \endtable

Usage
    \table
        Property      | Description                  | Req'd? | Default
        conductivity  | k single scalar value        | yes    |
        htc           | htc single scalar value      | yes    |
        Tinf          | Tinf single scalar value     | yes    |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        {
        type            convectionExternal;
        htc             3.7;
        conductivity    1;
        Tinf            4.2;
        value           uniform 0; // optional initial value
    }
    \endverbatim

SourceFiles
    adjntCOObjCOFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef adjntCOObjCOFvPatchScalarField_H
#define adjntCOObjCOFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class adjntCOObjCOFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class adjntCOObjCOFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data
    
        scalar TT_;

    // Private Member Functions

        //- Return current time
        scalar t() const;


public:

    //- Runtime type information
    TypeName("adjntCOObjCO");


    // Constructors

        //- Construct from patch and internal field
        adjntCOObjCOFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        adjntCOObjCOFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given mixedTypeFvPatchField
        //  onto a new patch
        adjntCOObjCOFvPatchScalarField
        (
            const adjntCOObjCOFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        adjntCOObjCOFvPatchScalarField
        (
            const adjntCOObjCOFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new adjntCOObjCOFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        adjntCOObjCOFvPatchScalarField
        (
            const adjntCOObjCOFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new adjntCOObjCOFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
