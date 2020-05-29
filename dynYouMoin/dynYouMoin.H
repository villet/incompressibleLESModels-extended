/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::incompressible::LESModels::dynYouMoin

Description
    You & Moin dynamic global Smagorinsky coefficient model.

SourceFiles
    dynYouMoin.C

\*---------------------------------------------------------------------------*/

#ifndef dynYouMoin_H
#define dynYouMoin_H

#include "Smagorinsky.H"
#include "LESfilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class YouMoin Declaration
\*---------------------------------------------------------------------------*/

class dynYouMoin
:
    public GenEddyVisc
{
    // Private data

        volScalarField k_;

        volScalarField Cs_;
        volScalarField PiSqr_;
        volScalarField PiFilteredSqr_;

        autoPtr<LESfilter> filterPtr_;
        LESfilter& filter_;


    // Private Member Functions

        //- Calculate coefficients cD, cI from filtering velocity field
        volScalarField PiSqr() const;
        volScalarField PiFilteredSqr() const;
        dimensionedScalar Cmu(const volSymmTensorField& D) const;
        dimensionedScalar cI(const volSymmTensorField& D) const;

        // Disallow default bitwise copy construct and assignment
        dynYouMoin(const dynYouMoin&);
        dynYouMoin& operator=(const dynYouMoin&);


public:

    //- Runtime type information
    TypeName("dynYouMoin");

    // Constructors

        //- Constructor from components
        dynYouMoin
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport
        );


    // Destructor

        ~dynYouMoin();


    // Member Functions

        //- Return SGS kinetic energy
        tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Correct Eddy-Viscosity and related properties
        void correct(const tmp<volTensorField>& gradU);

        //- Read turbulenceProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
