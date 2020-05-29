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

\*---------------------------------------------------------------------------*/

#include "dynOneEqEddy2.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dynOneEqEddy2, 0);
addToRunTimeSelectionTable(LESModel, dynOneEqEddy2, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

dimensionedScalar dynOneEqEddy2::ck(const volSymmTensorField& D) const
{
    volScalarField KK = 0.5*(filter_(magSqr(U())) - magSqr(filter_(U())));

    volSymmTensorField LL = dev(filter_(sqr(U())) - sqr(filter_(U())));

    volSymmTensorField MM =
        delta()*(filter_(sqrt(k_)*D) - 2*sqrt(KK + filter_(k_))*filter_(D));

    dimensionedScalar MMMM = average(magSqr(MM));

    if (MMMM.value() > VSMALL)
    {
        return average(LL && MM)/MMMM;
    }
    else
    {
        return 0.0;
    }
}


dimensionedScalar dynOneEqEddy2::ce(const volSymmTensorField& D) const
{
    volScalarField KK = 0.5*(filter_(magSqr(U())) - magSqr(filter_(U())));

    volScalarField mm =
        pow(KK + filter_(k_), 1.5)/(2*delta()) - filter_(pow(k_, 1.5))/delta();

    volScalarField ee =
        2*delta()*ck(D)
       *(
           filter_(sqrt(k_)*magSqr(D))
         - 2*sqrt(KK + filter_(k_))*magSqr(filter_(D))
        );

    dimensionedScalar mmmm = average(magSqr(mm));

    if (mmmm.value() > VSMALL)
    {
        return average(ee*mm)/mmmm;
    }
    else
    {
        return 0.0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dynOneEqEddy2::dynOneEqEddy2
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport
)
:
    LESModel(typeName, U, phi, transport),
    GenEddyVisc(U, phi, transport),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    Cs_
    (
        IOobject
        (
            "Cs",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    filterPtr_(LESfilter::New(U.mesh(), coeffDict())),
    filter_(filterPtr_())
{
    printCoeffs();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dynOneEqEddy2::~dynOneEqEddy2()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dynOneEqEddy2::correct(const tmp<volTensorField>& gradU)
{
    GenEddyVisc::correct(gradU);

    volSymmTensorField D = symm(gradU);

    volScalarField P = 2.0*nuSgs_*magSqr(D);

    solve
    (
       fvm::ddt(k_)
     + fvm::div(phi(), k_)
     - fvm::laplacian(DkEff(), k_)
    ==
       P
     - fvm::Sp(ce(D)*sqrt(k_)/delta(), k_)
    );

    bound(k_, k0());

    nuSgs_ = ck(D)*sqrt(k_)*delta();
    nuSgs_.correctBoundaryConditions();

    Cs_ = nuSgs_/(sqr(delta())*sqrt(magSqr(D)));

}


bool dynOneEqEddy2::read()
{
    if (GenEddyVisc::read())
    {
        filter_.read(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
