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

#include "dynSmagorinsky2.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dynSmagorinsky2, 0);
addToRunTimeSelectionTable(LESModel, dynSmagorinsky2, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

dimensionedScalar dynSmagorinsky2::cD(const volSymmTensorField& D) const
{
    volSymmTensorField LL = dev(filter_(sqr(U())) - (sqr(filter_(U()))));

    volSymmTensorField MM =
        sqr(delta())*(filter_(mag(D)*(D)) - 4*mag(filter_(D))*filter_(D));

    // kerrotaanko kahdella?
    dimensionedScalar MMMM = 2*average(magSqr(MM));

    if (MMMM.value() > VSMALL)
    {
        return average(LL && MM)/MMMM;
    }
    else
    {
        return 0.0;
    }
}


dimensionedScalar dynSmagorinsky2::cI(const volSymmTensorField& D) const
{
    volScalarField KK = 0.5*(filter_(magSqr(U())) - magSqr(filter_(U())));

    volScalarField mm =
        sqr(delta())*(4*sqr(mag(filter_(D))) - filter_(sqr(mag(D))));

    dimensionedScalar mmmm = average(magSqr(mm));

    if (mmmm.value() > VSMALL)
    {
        return average(KK*mm)/mmmm;
    }
    else
    {
        return 0.0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dynSmagorinsky2::dynSmagorinsky2
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
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("k",dimVelocity*dimVelocity, 0.0)
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

dynSmagorinsky2::~dynSmagorinsky2()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dynSmagorinsky2::correct(const tmp<volTensorField>& gradU)
{
    LESModel::correct(gradU);

    volSymmTensorField D = dev(symm(gradU));
    volScalarField magSqrD = magSqr(D);

    k_ = cI(D)*sqr(delta())*magSqrD;

    Cs_ = cD(D);

    nuSgs_ = Cs_*sqr(delta())*sqrt(magSqrD);
    nuSgs_.correctBoundaryConditions();
}


bool dynSmagorinsky2::read()
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
