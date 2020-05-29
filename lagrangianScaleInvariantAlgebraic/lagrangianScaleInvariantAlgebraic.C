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

#include "lagrangianScaleInvariantAlgebraic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(lagrangianScaleInvariantAlgebraic, 0);
addToRunTimeSelectionTable(LESModel, lagrangianScaleInvariantAlgebraic, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// !!! Should be updated every 5 time steps
volScalarField lagrangianScaleInvariantAlgebraic::cD(const volSymmTensorField& D)
{
Info<< "Entering function = " << runTime_.elapsedCpuTime() << " s." << endl;

  volSymmTensorField LL = dev(filter_(sqr(U())) - (sqr(filter_(U()))));

  volSymmTensorField MM =
    2*sqr(delta())*(filter_(mag(D)*(D)) - 4*mag(filter_(D))*filter_(D));

  if(initializeI_)
    {
      IMM_ = MM && MM;
      ILM_ = 0.16*0.16*IMM_; // initialization using standard smagorinsky coeff
      initializeI_ = false;
    }

  // !!! This should be moved to private data
  dimensionedScalar small
    (
     "small",
     dimensionSet(0, 8, -8, 0, 0, 0, 0),
     1e-9
   );

  // Relaxation time scale
  volScalarField T = theta_ * delta() * pow(ILM_ * IMM_ + small,-0.125);

  // Interpolation coefficient
  volScalarField epsilon = (runTime_.deltaT() / T) / (1 + runTime_.deltaT() / T);

  vectorField oldCoordinate = mesh_.C() - Uold_*runTime_.deltaT();
  Uold_ = U();

Info<< "Entering loop = " << runTime_.elapsedCpuTime() << " s." << endl;

  forAll(mesh_.cells(), celli)
  {
      label oldCell = mesh_.findCell(oldCoordinate[celli]);

      if (oldCell != -1)
      {
        // Calculate new integral values
  
       ILM_[celli] = epsilon[celli]*(LL[celli] && MM[celli]) 
              + (1-epsilon[celli])*ILM_[oldCell];
       IMM_[celli] = epsilon[celli]*(MM[celli] && MM[celli])
              + (1-epsilon[celli])*IMM_[oldCell];    
      }
  }

Info<< "Exit loop = " << runTime_.elapsedCpuTime() << " s." << endl;


   // !!! Should be moved to private data
   dimensionedScalar vsmall
     (
      "vsmall",
      dimensionSet(0, 4, -4, 0, 0, 0, 0),
      1e-32
      );
  
   // Values clipped to avoid complex values
   ILM_ = max(vsmall, ILM_);

   // Avoiding division by zero
   forAll(IMM_, celli)
     {
       if(IMM_[celli] < VSMALL)
	 {
	   IMM_[celli] = 1e-9;
	 }
     }

Info<< "Exit function = " << runTime_.elapsedCpuTime() << " s." << endl;

   return (ILM_ / IMM_);
}

  // Assumed same global dynamic procedure as in "dynSmagorinky" model
  dimensionedScalar lagrangianScaleInvariantAlgebraic::cI(const volSymmTensorField& D) const
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

lagrangianScaleInvariantAlgebraic::lagrangianScaleInvariantAlgebraic
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

    Uold_(U),

    filterPtr_(LESfilter::New(U.mesh(), coeffDict())),
    filter_(filterPtr_())
{
    printCoeffs();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

lagrangianScaleInvariantAlgebraic::~lagrangianScaleInvariantAlgebraic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void lagrangianScaleInvariantAlgebraic::correct(const tmp<volTensorField>& gradU)
{
  LESModel::correct(gradU);

  volSymmTensorField D = dev(symm(gradU));
  volScalarField magSqrD = magSqr(D);

  k_ = cI(D)*sqr(delta())*magSqrD;

  nuSgs_ = cD(D)*sqr(delta())*sqrt(magSqrD);
  nuSgs_.correctBoundaryConditions();
}


bool lagrangianScaleInvariantAlgebraic::read()
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
