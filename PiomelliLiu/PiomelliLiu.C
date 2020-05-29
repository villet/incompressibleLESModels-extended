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

// USE BOUND FUNCTION

#include "PiomelliLiu.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(PiomelliLiu, 0);
addToRunTimeSelectionTable(LESModel, PiomelliLiu, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

volScalarField PiomelliLiu::cD(const volSymmTensorField& D)
{

  volSymmTensorField LL = dev(filter_(sqr(U())) - (sqr(filter_(U()))));
  volSymmTensorField alpha = 4*filter_(sqr(delta()))*mag(filter_(D))*filter_(D);
  volSymmTensorField beta = sqr(delta())*mag(D)*D;

  volScalarField denominator = alpha && alpha;
  volScalarField numerator = (LL - 2*filter_((Cs_+runTime_.deltaT()/deltaTold_*(CsOld_- 
                             CsOldOld_))*beta)) && alpha;

  CsOldOld_ = CsOld_;
  CsOld_ = Cs_;

  deltaTold_ = runTime_.deltaT();
 
  return (-0.5*numerator/denominator);
}

  // Used global ksgs coefficient "dynSmagorinky" model
  dimensionedScalar PiomelliLiu::cI(const volSymmTensorField& D) const
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
  

void PiomelliLiu::clipMin(const volScalarField& ref, volScalarField& corr)
{
  label clipped = 0;

  forAll(corr, celli)
  {
    if (corr[celli] < ref[celli])
      {
	corr[celli] = ref[celli] + SMALL;
	clipMin_[celli] = clipMin_[celli] + 1.0;
	clipped++;
      }
  }

  Info << "Clipped " << clipped << " too low Cs values." << endl;
}

void PiomelliLiu::clipMax(const scalar& ref, volScalarField& corr)
{
  label clipped = 0;

  forAll(corr, celli)
  {
    if (corr[celli] > ref)
    {
      corr[celli] = ref;
      clipMax_[celli] = clipMin_[celli] + 1.0;
      clipped++;
    }
  }

  Info << "Clipped " << clipped << " Cs values higher than 1.0^2." << endl; 
 
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

PiomelliLiu::PiomelliLiu
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

    clipMax_
    (
        IOobject
        (
            "clipMax",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
	dimensionedScalar("clipMax", dimless, 0.0)
    ),

    clipMin_
    (
        IOobject
        (
            "clipMin",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
	dimensionedScalar("clipMin", dimless, 0.0)
    ),

    allowBackscatter_
    (
        Switch::lookupOrAddToDict
        (
            "allowBackscatter",
            coeffDict(),
            false
        )
    ),

    CsOld_ (Cs_),
    CsOldOld_ (CsOld_),

    deltaTold_ (runTime_.deltaT()),

    filterPtr_(LESfilter::New(U.mesh(), coeffDict())),
    filter_(filterPtr_())
{
    printCoeffs();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

PiomelliLiu::~PiomelliLiu()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void PiomelliLiu::correct(const tmp<volTensorField>& gradU)
{
  LESModel::correct(gradU);

  volSymmTensorField D = dev(symm(gradU));
  volScalarField magSqrD = magSqr(D);

  k_ = cI(D)*sqr(delta())*magSqrD;

  Cs_ = cD(D);

  if (allowBackscatter_)
  {
    clipMin(-nu()/(sqr(delta())*sqrt(magSqrD)), Cs_);
    clipMax(0.2*0.2, Cs_);
  }
  else
  {
    clipMin(Cs_*0.0, Cs_);
    clipMax(1.0*1.0, Cs_);
  }
  Cs_.correctBoundaryConditions();

  nuSgs_ = Cs_*sqr(delta())*sqrt(magSqrD);
  nuSgs_.correctBoundaryConditions();

}


bool PiomelliLiu::read()
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
