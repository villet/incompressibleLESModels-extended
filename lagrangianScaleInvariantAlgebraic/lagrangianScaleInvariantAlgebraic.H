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
    Foam::incompressible::LESModels::lagrangianScaleInvariantAlgebraic

Description
    Lagrangian-averaged dynamic Smagorinsky Model for incompressible flows.

SourceFiles
    lagrangianScaleInvariantAlgebraic.C

\*---------------------------------------------------------------------------*/

#ifndef lagrangianScaleInvariantAlgebraic_H
#define lagrangianScaleInvariantAlgebraic_H

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
                           Class lagrangianScaleInvariant Declaration
\*---------------------------------------------------------------------------*/

class lagrangianScaleInvariantAlgebraic
:
    public GenEddyVisc
{
    // Private data

        volScalarField k_;
        volScalarField ILM_;
        volScalarField IMM_;
        dimensionedScalar theta_;
        Switch initializeI_;
        volVectorField Uold_;

        autoPtr<LESfilter> filterPtr_;
        LESfilter& filter_;


    // Private Member Functions

        //- Calculate coefficients cD, cI from filtering velocity field
        volScalarField cD(const volSymmTensorField& D);
        dimensionedScalar cI(const volSymmTensorField& D) const;

        // Disallow default bitwise copy construct and assignment
        lagrangianScaleInvariantAlgebraic(const lagrangianScaleInvariantAlgebraic&);
        lagrangianScaleInvariantAlgebraic& operator=(const lagrangianScaleInvariantAlgebraic&);


public:

    //- Runtime type information
    TypeName("lagrangianScaleInvariantAlgebraic");

    // Constructors

        //- Constructor from components
        lagrangianScaleInvariantAlgebraic
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport
        );


    // Destructor

        ~lagrangianScaleInvariantAlgebraic();


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
