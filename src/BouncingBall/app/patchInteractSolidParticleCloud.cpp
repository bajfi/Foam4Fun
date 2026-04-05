// Copyright (c) 2026 JackLee
//
// This software is released under the MIT License.
// https://opensource.org/licenses/MIT

#include "patchInteractSolidParticleCloud.h"
#include "fvMesh.H"
#include "volFields.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchInteractSolidParticleCloud::patchInteractSolidParticleCloud(
  const fvMesh& mesh, const word& cloudName, bool readFields
)
  : Cloud<patchInteractSolidParticle>(mesh, cloudName, false),
    mesh_(mesh),
    particleProperties_(IOobject(
      "particleProperties",
      mesh_.time().constant(),
      mesh_,
      IOobject::MUST_READ_IF_MODIFIED,
      IOobject::NO_WRITE
    )),
    physicalProperties_{
      dimensionedScalar("rhop", particleProperties_).value(),
      dimensionedScalar("e", particleProperties_).value(),
      dimensionedScalar("mu", particleProperties_).value()
    }
{
    if (physicalProperties_.e < 0 || physicalProperties_.e > 1)
    {
        FatalErrorInFunction << "Coefficient of restitution e must be between 0 and 1, but got "
                             << physicalProperties_.e << exit(FatalError);
    }

    if (readFields)
    {
        patchInteractSolidParticle::readFields(*this);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::patchInteractSolidParticleCloud::move(const dimensionedVector& g)
{
    const volScalarField& rho = mesh_.lookupObject<const volScalarField>("rho");
    const volVectorField& U = mesh_.lookupObject<const volVectorField>("U");
    const volScalarField& nu = mesh_.lookupObject<const volScalarField>("nu");

    interpolationCellPoint<scalar> rhoInterp(rho);
    interpolationCellPoint<vector> UInterp(U);
    interpolationCellPoint<scalar> nuInterp(nu);

    patchInteractSolidParticle::trackingData td(
      *this, rhoInterp, UInterp, nuInterp, g.value(), physicalProperties_
    );

    Cloud<patchInteractSolidParticle>::move(*this, td, mesh_.time().deltaTValue());
}

// ************************************************************************* //
