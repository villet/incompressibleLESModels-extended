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
    Foam::incompressible::LESModels::Smagorinsky3

Description
    The Isochoric Smagorinsky3 Model for incompressible flows.

    Algebraic eddy viscosity SGS model founded on the assumption that
    local equilibrium prevails.
    Thus,
    @verbatim
        B = 2/3*k*I - 2*nuSgs*dev(D)
        Beff = 2/3*k*I - 2*nuEff*dev(D)

    where

        D = symm(grad(U));
        k = (2*ck/ce)*delta^2*||D||^2
        nuSgs = ck*sqrt(k)*delta
        nuEff = nuSgs + nu
    @endverbatim

SourceFiles
    Smagorinsky3.C

\*---------------------------------------------------------------------------*/

#ifndef Smagorinsky3_H
#define Smagorinsky3_H

#include "GenEddyVisc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class Smagorinsky3 Declaration
\*---------------------------------------------------------------------------*/

class Smagorinsky3
:
    public GenEddyVisc
{
    // Private data

        dimensionedScalar ck_;
        dimensionedScalar Cs_;

    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        Smagorinsky3(const Smagorinsky3&);
        Smagorinsky3& operator=(const Smagorinsky3&);


public:

    //- Runtime type information
    TypeName("Smagorinsky3");

    // Constructors

        //- Construct from components
        Smagorinsky3
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport
        );


    // Destructor

        ~Smagorinsky3()
        {}


    // Member Functions

        //- Return SGS kinetic energy
        //  calculated from the given velocity gradient
        tmp<volScalarField> k(const tmp<volTensorField>& gradU) const
        {
            return (2.0*ck_/ce_)*sqr(delta())*magSqr(dev(symm(gradU)));
        }

        //- Return SGS kinetic energy
        tmp<volScalarField> k() const
        {
            return k(fvc::grad(U()));
        }


        //- Correct Eddy-Viscosity and related properties
        virtual void correct(const tmp<volTensorField>& gradU);

        //- Read turbulenceProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
