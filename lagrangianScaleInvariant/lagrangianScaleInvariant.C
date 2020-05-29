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

// DO NOT REQUIRE ILM & IMM FIELDS (NO WRITE)
// USE BOUND FUNCTION
// CALCULATE KSGS USING LOCALIZED MODEL

#include "lagrangianScaleInvariant.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(lagrangianScaleInvariant, 0);
addToRunTimeSelectionTable(LESModel, lagrangianScaleInvariant, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

volScalarField lagrangianScaleInvariant::cD(const volSymmTensorField& D)
{
  volSymmTensorField LL = dev(filter_(sqr(U())) - (sqr(filter_(U()))));

  volSymmTensorField MM =
    sqr(delta())*(filter_(mag(D)*(D)) - 4*mag(filter_(D))*filter_(D));

  if(initializeI_)
    {
      IMM_ = MM && MM;
      ILM_ = 0.16*0.16*IMM_; // initialization using standard Smagorinsky coeff
      initializeI_ = false;
    }

  // This should be moved to private data
  dimensionedScalar small
    (
     "small",
     dimensionSet(0, 8, -8, 0, 0, 0, 0),
     1e-9
   );

  // Relaxation time scale
  volScalarField T = theta_ * delta() * pow(ILM_ * IMM_ + small,-0.125);

  solve
  (
    fvm::ddt(ILM_)
    + fvm::div(phi(),ILM_)
    ==
    (LL && MM)/T 
    - fvm::Sp(1/T,ILM_)
  );

  solve
  (
    fvm::ddt(IMM_)
    + fvm::div(phi(),IMM_)
    ==
    (MM && MM)/T 
    - fvm::Sp(1/T,IMM_)
  );
      
  dimensionedScalar vsmall
    (
     "vsmall",
     dimensionSet(0, 4, -4, 0, 0, 0, 0),
     1e-32
     );
   
   // Values clipped to avoid complex values
   ILM_ = max(vsmall, ILM_);

   // Avoiding division by zero, not really 
   forAll(IMM_, celli)
     {
       if(IMM_[celli] < VSMALL)
	 {
	   IMM_[celli] = VSMALL;
	 }
     }

   return (ILM_ / IMM_);
}

  // Assumed same global ksgs dynamic procedure as in "dynSmagorinky" model
  dimensionedScalar lagrangianScaleInvariant::cI(const volSymmTensorField& D) const
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

lagrangianScaleInvariant::lagrangianScaleInvariant
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

    ILM_
    (
        IOobject
        (
            "ILM",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    IMM_
    (
        IOobject
        (
            "IMM",
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

    theta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "theta",
            coeffDict(),
            1.5
        )
    ),

    initializeI_
    (
        Switch::lookupOrAddToDict
        (
            "initializeI",
            coeffDict(),
            false
        )
    ),

    filterPtr_(LESfilter::New(U.mesh(), coeffDict())),
    filter_(filterPtr_())
{
    printCoeffs();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

lagrangianScaleInvariant::~lagrangianScaleInvariant()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void lagrangianScaleInvariant::correct(const tmp<volTensorField>& gradU)
{
  LESModel::correct(gradU);

  volSymmTensorField D = dev(symm(gradU));
  volScalarField magSqrD = magSqr(D);

  k_ = cI(D)*sqr(delta())*magSqrD;

  Cs_ = cD(D);

  nuSgs_ = Cs_*sqr(delta())*sqrt(magSqrD);
  nuSgs_.correctBoundaryConditions();
}


bool lagrangianScaleInvariant::read()
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
