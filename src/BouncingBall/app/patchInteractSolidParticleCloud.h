// Copyright (c) 2026 JackLee
//
// This software is released under the MIT License.
// https://opensource.org/licenses/MIT

#pragma once

#include "Cloud.H"
#include "patchInteractSolidParticle.h"
#include "dimensionedVector.H"

namespace Foam
{

class fvMesh;

class patchInteractSolidParticleCloud : public Cloud<patchInteractSolidParticle>
{
    // Private Data

    const fvMesh& mesh_;

    IOdictionary particleProperties_;

    patchInteractPhysicalProperties physicalProperties_;

  public:
    //- Type of parcel within the cloud
    typedef patchInteractSolidParticle parcelType;

    //- No copy construct
    patchInteractSolidParticleCloud(const patchInteractSolidParticleCloud&) = delete;

    //- No copy assignment
    void operator=(const patchInteractSolidParticleCloud&) = delete;

    // Constructors

    //- Read construct
    explicit patchInteractSolidParticleCloud(
      const fvMesh& mesh, const word& cloudName = cloud::defaultName, bool readFields = true
    );

    // Member Functions

    // Access

    const fvMesh& mesh() const
    {
        return mesh_;
    }

    const patchInteractPhysicalProperties& physicalProperties() const noexcept
    {
        return physicalProperties_;
    }

    // Edit

    //- Move the particles under the influence of the given
    //-  gravitational acceleration
    void move(const dimensionedVector& g);
};
} // namespace Foam
