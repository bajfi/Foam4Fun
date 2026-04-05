// Copyright (c) 2026 JackLee
//
// This software is released under the MIT License.
// https://opensource.org/licenses/MIT
#pragma once

#include "particle.H"
#include "interpolationCellPoint.H"
#include "PatchInteractionModel.H"

namespace Foam
{

class patchInteractSolidParticle;

struct patchInteractPhysicalProperties
{
    scalar rhop;
    scalar e;
    scalar mu;
};

// Namespace-scope declaration so operator<< is declared within Foam
Ostream& operator<<(Ostream&, const patchInteractSolidParticle&);

class patchInteractSolidParticle : public particle
{
    // Private Data

    //- Diameter
    scalar d_;

    //- Velocity of parcel
    vector U_;

  public:
    friend class Cloud<patchInteractSolidParticle>;

    //- Class used to pass tracking data to the trackToFace function
    class trackingData : public particle::trackingData
    {
        // Interpolators for continuous phase fields

        const interpolationCellPoint<scalar>& rhoInterp_;
        const interpolationCellPoint<vector>& UInterp_;
        const interpolationCellPoint<scalar>& nuInterp_;

        //- Local gravitational or other body-force acceleration
        const vector& g_;

        //- Particle interaction/transport physical properties
        const patchInteractPhysicalProperties& physicalProperties_;

      public:
        // Constructors

        trackingData(
          const Cloud<patchInteractSolidParticle>& cloud,
          const interpolationCellPoint<scalar>& rhoInterp,
          const interpolationCellPoint<vector>& UInterp,
          const interpolationCellPoint<scalar>& nuInterp,
          const vector& g,
          const patchInteractPhysicalProperties& physicalProperties
        )
          : particle::trackingData(cloud),
            rhoInterp_(rhoInterp),
            UInterp_(UInterp),
            nuInterp_(nuInterp),
            g_(g),
            physicalProperties_(physicalProperties)
        {}

        // Member Functions

        const interpolationCellPoint<scalar>& nuInterp() const noexcept
        {
            return nuInterp_;
        }

        const interpolationCellPoint<vector>& UInterp() const noexcept
        {
            return UInterp_;
        }

        const interpolationCellPoint<scalar>& rhoInterp() const noexcept
        {
            return rhoInterp_;
        }

        const vector& g() const noexcept
        {
            return g_;
        }

        const patchInteractPhysicalProperties& physicalProperties() const noexcept
        {
            return physicalProperties_;
        }
    };

    // Static data members

    //- Size in bytes of the fields
    static const std::size_t sizeofFields;

    // Constructors

    //- Construct from a position and a cell
    //  Searches for the rest of the required topology.
    //  Other properties are zero initialised.
    patchInteractSolidParticle(const polyMesh& mesh, const vector& position, const label celli = -1)
      : particle(mesh, position, celli), d_(0), U_(Zero)
    {}

    //- Construct from components
    patchInteractSolidParticle(
      const polyMesh& mesh,
      const barycentric& coordinates,
      const label celli,
      const label tetFacei,
      const label tetPti,
      const scalar d,
      const vector& U
    )
      : particle(mesh, coordinates, celli, tetFacei, tetPti), d_(d), U_(U)
    {}

    //- Construct from Istream
    patchInteractSolidParticle(
      const polyMesh& mesh, Istream& is, bool readFields = true, bool newFormat = true
    );

    //- Construct and return a clone
    virtual autoPtr<particle> clone() const
    {
        return autoPtr<particle>(new patchInteractSolidParticle(*this));
    }

    //- Factory class to read-construct particles (for parallel transfer)
    class iNew
    {
        const polyMesh& mesh_;

      public:
        iNew(const polyMesh& mesh) : mesh_(mesh) {}

        autoPtr<patchInteractSolidParticle> operator()(Istream& is) const
        {
            return autoPtr<patchInteractSolidParticle>::New(mesh_, is, true);
        }
    };

    // Member Functions

    // Access

    //- Return diameter
    scalar d() const noexcept
    {
        return d_;
    }

    //- Return velocity
    const vector& U() const noexcept
    {
        return U_;
    }

    // Tracking

    //- Move
    bool move(Cloud<patchInteractSolidParticle>&, trackingData&, const scalar);

    // Patch interactions

    //- Overridable function to handle the particle hitting a patch
    //  Executed before other patch-hitting functions
    bool hitPatch(Cloud<patchInteractSolidParticle>& cloud, trackingData& td);

    //- Overridable function to handle the particle hitting a
    //  processorPatch
    void hitProcessorPatch(Cloud<patchInteractSolidParticle>& cloud, trackingData& td);

    //- Overridable function to handle the particle hitting a wallPatch
    void hitWallPatch(Cloud<patchInteractSolidParticle>& cloud, trackingData& td);

    // Physics core (cloud-independent)

    //- Pure kernel: compute updated particle velocity from drag+body force
    static vector computeVelocityAfterDrag(
      const vector& Up,
      const scalar d,
      const scalar dt,
      const scalar rhoc,
      const vector& Uc,
      const scalar nuc,
      const vector& g,
      const patchInteractPhysicalProperties& physicalProperties
    );

    //- Pure kernel: compute post-wall-interaction velocity
    static vector computeVelocityAfterWallInteraction(
      const vector& Up,
      const vector& wallNormal,
      const patchInteractPhysicalProperties& physicalProperties
    );

    //- Advance particle velocity over dt using local carrier properties
    void updateVelocity(
      const scalar dt,
      const scalar rhoc,
      const vector& Uc,
      const scalar nuc,
      const vector& g,
      const patchInteractPhysicalProperties& physicalProperties
    );

    //- Apply wall interaction response using a supplied wall normal
    void applyWallInteraction(
      const vector& wallNormal, const patchInteractPhysicalProperties& physicalProperties
    );

    //- Transform the physical properties of the particle
    //  according to the given transformation tensor
    virtual void transformProperties(const tensor& T);

    //- Transform the physical properties of the particle
    //  according to the given separation vector
    virtual void transformProperties(const vector& separation);

    // I-O

    static void readFields(Cloud<patchInteractSolidParticle>& c);

    static void writeFields(const Cloud<patchInteractSolidParticle>& c);

    // Ostream Operator

    friend Ostream& operator<<(Ostream&, const patchInteractSolidParticle&);
};

// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Contiguous data for patchInteractSolidParticle
template <>
struct is_contiguous<patchInteractSolidParticle> : std::true_type
{};

} // namespace Foam
