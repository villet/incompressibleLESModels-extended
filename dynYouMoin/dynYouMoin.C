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

#include "dynYouMoin.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dynYouMoin, 0);
addToRunTimeSelectionTable(LESModel, dynYouMoin, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

volScalarField dynYouMoin::PiSqr() const
{
  volTensorField alpha = fvc::grad(U());

  volVectorField alphaVecX = vector(1,0,0) & alpha;
  volVectorField alphaVecY = vector(0,1,0) & alpha;
  volVectorField alphaVecZ = vector(0,0,1) & alpha;

  volTensorField beta = sqr(delta())*(alphaVecX * alphaVecX + alphaVecY * alphaVecY + alphaVecZ * alphaVecZ);

//  volTensorField beta = sqr(delta())*alpha & alpha;
  volScalarField B = beta.component(tensor::XX)*beta.component(tensor::YY)
                   - beta.component(tensor::XY)*beta.component(tensor::XY)
                   + beta.component(tensor::XX)*beta.component(tensor::ZZ) 
                   - beta.component(tensor::XZ)*beta.component(tensor::XZ)
                   + beta.component(tensor::YY)*beta.component(tensor::ZZ)
                   - beta.component(tensor::YZ)*beta.component(tensor::YZ);

return B/(alpha && alpha);
}

volScalarField dynYouMoin::PiFilteredSqr() const
{
  volTensorField alpha = fvc::grad(filter_(U()));

  volVectorField alphaVecX = vector(1,0,0) & alpha;
  volVectorField alphaVecY = vector(0,1,0) & alpha;
  volVectorField alphaVecZ = vector(0,0,1) & alpha;

  volTensorField beta = 4*sqr(delta())*(alphaVecX * alphaVecX + alphaVecY * alphaVecY + alphaVecZ * alphaVecZ);

//  volTensorField beta = 4*sqr(filter_(delta()))*alpha & alpha;
  volScalarField B = beta.component(tensor::XX)*beta.component(tensor::YY)
                   - beta.component(tensor::XY)*beta.component(tensor::XY)
                   + beta.component(tensor::XX)*beta.component(tensor::ZZ) 
                   - beta.component(tensor::XZ)*beta.component(tensor::XZ)
                   + beta.component(tensor::YY)*beta.component(tensor::ZZ)
                   - beta.component(tensor::YZ)*beta.component(tensor::YZ);

  return B / (alpha && alpha);
}

dimensionedScalar dynYouMoin::Cmu(const volSymmTensorField& D) const
{
  volTensorField alphaFiltered = fvc::grad(filter_(U()));
  volTensorField alpha = fvc::grad(U());

  dimensionedScalar numerator = average(filter_(magSqr(alpha))
                              - magSqr(alphaFiltered));

  dimensionedScalar denominator = average(filter_(sqrt(PiSqr_) * magSqr(D))
                                - sqrt(PiFilteredSqr_) * magSqr(filter_(D)));

  dimensionedScalar result = -0.5*average(nu())*numerator/denominator;

  if (result.value() > VSMALL)
    {
      return result;
    }
  else
    {
      return 0.0;
    }
}

  // Used global ksgs coefficient "dynSmagorinky" model
  dimensionedScalar dynYouMoin::cI(const volSymmTensorField& D) const
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

dynYouMoin::dynYouMoin
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

    PiSqr_
    (
        IOobject
        (
            "PiSqr",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
	dimensionedScalar("PiSqr",dimensionSet(0,4,-2,0,0,0,0), 0.0)
    ),

    PiFilteredSqr_
    (
        IOobject
        (
            "PiFilteredSqr",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
	dimensionedScalar("PiFilteredSqr",dimensionSet(0,4,-2,0,0,0,0), 0.0)
    ),

    filterPtr_(LESfilter::New(U.mesh(), coeffDict())),
    filter_(filterPtr_())
{
    printCoeffs();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dynYouMoin::~dynYouMoin()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dynYouMoin::correct(const tmp<volTensorField>& gradU)
{
  LESModel::correct(gradU);

  volSymmTensorField D = dev(symm(gradU));
  volScalarField magSqrD = magSqr(D);

  k_ = cI(D)*sqr(delta())*magSqrD;

  PiSqr_ = PiSqr();
  bound(PiSqr_, dimensionedScalar("PiSqrBound",dimensionSet(0,4,-2,0,0,0,0), VSMALL));

  PiFilteredSqr_ = PiFilteredSqr();
  bound(PiFilteredSqr_, dimensionedScalar("PiFilteredSqrBound",dimensionSet(0,4,-2,0,0,0,0), VSMALL));

  nuSgs_ = Cmu(D)*sqrt(PiSqr_);
  nuSgs_.correctBoundaryConditions();

  Cs_ = nuSgs_/(sqr(delta())*sqrt(magSqr(D)));

}


bool dynYouMoin::read()
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
