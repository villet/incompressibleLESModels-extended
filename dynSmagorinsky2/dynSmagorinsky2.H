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
    Foam::incompressible::LESModels::dynSmagorinsky2

Description
    The Isochoric dynamic Smagorinsky Model for incompressible flows.

    Algebraic eddy viscosity SGS model founded on the assumption that
    local equilibrium prevails.
    Thus,
    @verbatim
        B = 2/3*k*I - 2*nuSgs*dev(D)
        Beff = 2/3*k*I - 2*nuEff*dev(D)

    where

        k = cI*delta^2*||D||^2
        nuSgs = cD*delta^2*||D||
        nuEff = nuSgs + nu

    In the dynamic version of the choric  Smagorinsky model
    the coefficients cI and cD are calculated during the simulation,

        cI=<K*m>/<m*m>

    and

        cD=<L.M>/<M.M>,

    where

        K = 0.5*(F(U.U) - F(U).F(U))
        m = delta^2*(4*||F(D)||^2 - F(||D||^2))
        L = dev(F(U*U) - F(U)*F(U))
        M = delta^2*(F(||D||*dev(D)) - 4*||F(D)||*F(dev(D)))
    @endverbatim

SourceFiles
    dynSmagorinsky2.C

\*---------------------------------------------------------------------------*/

#ifndef dynSmagorinsky2_H
#define dynSmagorinsky2_H

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
                           Class dynSmagorinsky2 Declaration
\*---------------------------------------------------------------------------*/

class dynSmagorinsky2
:
    public GenEddyVisc
{
    // Private data

        volScalarField k_;
        volScalarField Cs_;

        autoPtr<LESfilter> filterPtr_;
        LESfilter& filter_;


    // Private Member Functions

        //- Calculate coefficients cD, cI from filtering velocity field
        dimensionedScalar cD(const volSymmTensorField& D) const;
        dimensionedScalar cI(const volSymmTensorField& D) const;

        // Disallow default bitwise copy construct and assignment
        dynSmagorinsky2(const dynSmagorinsky2&);
        dynSmagorinsky2& operator=(const dynSmagorinsky2&);


public:

    //- Runtime type information
    TypeName("dynSmagorinsky2");

    // Constructors

        //- Constructor from components
        dynSmagorinsky2
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport
        );


    // Destructor

        ~dynSmagorinsky2();


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
