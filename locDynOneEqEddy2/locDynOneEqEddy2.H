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
    Foam::incompressible::LESModels::locDynOneEqEddy2

Description
    Localised Dynamic One Equation Eddy Viscosity Model for incompressible
    flows

    Eddy viscosity SGS model using a modeled balance equation to simulate
    the behaviour of k, hence
    @verbatim
        d/dt(k) + div(U*k) - div(nuSgs*grad(k))
        =
        -B*L - ce*rho*k^3/2/delta
    and
        B = 2/3*k*I - 2*nuSgs*dev(D)
        Beff = 2/3*k*I - 2*nuEff*dev(D)
    where
        nuSgs = cD*delta^2*||D||
        nuEff = nuSgs + nu
    @endverbatim

    A dynamic procedure is here applied to evaluate ck and ce
    @verbatim
        ck=<L.M>/<M.M>
    and
        ce=<e*m>/<m*m>
    where
       K = 0.5*(F(U.U) - F(U).F(U))
       L = (F(U*U) - F(U)*F(U) - 0.33*K*I)
       M = delta*(F(sqrt(k)*D) - 2*sqrt(K + filter(k))*F(D))
       m = pow(K + F(k), 3.0/2.0)/(2*delta) - F(pow(k, 3.0/2.0))/delta
       e = 2*delta*ck*(F(sqrt(k)*(D && D)) - 2*sqrt(K + F(k))*(F(D) && F(D)))/
    @endverbatim

SourceFiles
    locDynOneEqEddy2.C

\*---------------------------------------------------------------------------*/

#ifndef locDynOneEqEddy2_H
#define locDynOneEqEddy2_H

#include "GenEddyVisc.H"
#include "simpleFilter.H"
#include "LESfilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class locDynOneEqEddy2 Declaration
\*---------------------------------------------------------------------------*/

class locDynOneEqEddy2
:
    public GenEddyVisc
{
    // Private data

        volScalarField k_;

        volScalarField Cs_;

        simpleFilter simpleFilter_;
        autoPtr<LESfilter> filterPtr_;
        LESfilter& filter_;


    // Private Member Functions

        //- Calculate ck, ce by filtering the velocity field U.
        volScalarField ck
        (
            const volSymmTensorField&,
            const volScalarField&
        ) const;

        volScalarField ce
        (
            const volSymmTensorField&,
            const volScalarField&
        ) const;

        // Disallow default bitwise copy construct and assignment
        locDynOneEqEddy2(const locDynOneEqEddy2&);
        locDynOneEqEddy2& operator=(const locDynOneEqEddy2&);


public:

    //- Runtime type information
    TypeName("locDynOneEqEddy2");

    // Constructors

        //- from components
        locDynOneEqEddy2
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport
        );


    // Destructor

        ~locDynOneEqEddy2();


    // Member Functions

        //- Return SGS kinetic energy
        tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", nuSgs_ + nu())
            );
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
