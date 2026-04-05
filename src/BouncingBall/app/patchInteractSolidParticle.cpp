// Copyright (c) 2026 JackLee
//
// This software is released under the MIT License.
// https://opensource.org/licenses/MIT

#include "patchInteractSolidParticle.h"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
defineTemplateTypeNameAndDebug(Cloud<patchInteractSolidParticle>, 0);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::patchInteractSolidParticle::move(
  Cloud<patchInteractSolidParticle>& cloud, trackingData& td, const scalar trackTime
)
{
    td.keepParticle = true;
    td.switchProcessor = false;

    while (td.keepParticle && !td.switchProcessor && stepFraction() < 1)
    {
        const scalar sfrac = stepFraction();
        const scalar f = 1 - stepFraction();

        trackToAndHitFace(f * trackTime * U_, f, cloud, td);

        const scalar dt = (stepFraction() - sfrac) * trackTime;

        const tetIndices tetIs = this->currentTetIndices();
        scalar rhoc = td.rhoInterp().interpolate(this->coordinates(), tetIs);
        vector Uc = td.UInterp().interpolate(this->coordinates(), tetIs);
        scalar nuc = td.nuInterp().interpolate(this->coordinates(), tetIs);

        updateVelocity(dt, rhoc, Uc, nuc, td.g(), td.physicalProperties());
    }

    return td.keepParticle;
}

bool Foam::patchInteractSolidParticle::hitPatch(Cloud<patchInteractSolidParticle>&, trackingData&)
{
    return false;
}

void Foam::patchInteractSolidParticle::hitProcessorPatch(
  Cloud<patchInteractSolidParticle>&, trackingData& td
)
{
    td.switchProcessor = true;
}

void Foam::patchInteractSolidParticle::hitWallPatch(
  Cloud<patchInteractSolidParticle>&, trackingData& td
)
{
    applyWallInteraction(normal(), td.physicalProperties());
}

void Foam::patchInteractSolidParticle::updateVelocity(
  const scalar dt,
  const scalar rhoc,
  const vector& Uc,
  const scalar nuc,
  const vector& g,
  const patchInteractPhysicalProperties& physicalProperties
)
{
    U_ = computeVelocityAfterDrag(U_, d_, dt, rhoc, Uc, nuc, g, physicalProperties);
}

void Foam::patchInteractSolidParticle::applyWallInteraction(
  const vector& wallNormal, const patchInteractPhysicalProperties& physicalProperties
)
{
    U_ = computeVelocityAfterWallInteraction(U_, wallNormal, physicalProperties);
}

Foam::vector Foam::patchInteractSolidParticle::computeVelocityAfterDrag(
  const vector& Up,
  const scalar d,
  const scalar dt,
  const scalar rhoc,
  const vector& Uc,
  const scalar nuc,
  const vector& g,
  const patchInteractPhysicalProperties& physicalProperties
)
{
    const scalar rhop = physicalProperties.rhop;
    const scalar magUr = mag(Uc - Up);

    scalar ReFunc = 1.0;
    const scalar Re = magUr * d / nuc;

    if (Re > 0.01)
    {
        ReFunc += 0.15 * pow(Re, 0.687);
    }

    const scalar Dc = (24.0 * nuc / d) * ReFunc * (3.0 / 4.0) * (rhoc / (d * rhop));

    return (Up + dt * (Dc * Uc + (1.0 - rhoc / rhop) * g)) / (1.0 + dt * Dc);
}

Foam::vector Foam::patchInteractSolidParticle::computeVelocityAfterWallInteraction(
  const vector& Up,
  const vector& wallNormal,
  const patchInteractPhysicalProperties& physicalProperties
)
{
    vector U = Up;

    const scalar Un = U & wallNormal;
    const vector Ut = U - Un * wallNormal;

    if (Un > 0)
    {
        U -= (1.0 + physicalProperties.e) * Un * wallNormal;
    }

    U -= physicalProperties.mu * Ut;
    return U;
}

void Foam::patchInteractSolidParticle::transformProperties(const tensor& T)
{
    particle::transformProperties(T);
    U_ = transform(T, U_);
}

void Foam::patchInteractSolidParticle::transformProperties(const vector& separation)
{
    particle::transformProperties(separation);
}

const std::size_t Foam::patchInteractSolidParticle::sizeofFields(
  sizeof(patchInteractSolidParticle) - sizeof(particle)
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchInteractSolidParticle::patchInteractSolidParticle(
  const polyMesh& mesh, Istream& is, bool readFields, bool newFormat
)
  : particle(mesh, is, readFields, newFormat)
{
    if (readFields)
    {
        if (is.format() == IOstreamOption::ASCII)
        {
            is >> d_ >> U_;
        }
        else if (!is.checkLabelSize<>() || !is.checkScalarSize<>())
        {
            // Non-native label or scalar size

            is.beginRawRead();

            readRawScalar(is, &d_);
            readRawScalar(is, U_.data(), vector::nComponents);

            is.endRawRead();
        }
        else
        {
            is.read(reinterpret_cast<char*>(&d_), sizeofFields);
        }
    }

    is.check(FUNCTION_NAME);
}

void Foam::patchInteractSolidParticle::readFields(Cloud<patchInteractSolidParticle>& c)
{
    const bool readOnProc = c.size();

    particle::readFields(c);

    IOField<scalar> d(c.fieldIOobject("d", IOobject::MUST_READ), readOnProc);
    c.checkFieldIOobject(c, d);

    IOField<vector> U(c.fieldIOobject("U", IOobject::MUST_READ), readOnProc);
    c.checkFieldIOobject(c, U);

    label i = 0;
    for (auto& p : c)
    {
        p.d_ = d[i];
        p.U_ = U[i];
        ++i;
    }
}

void Foam::patchInteractSolidParticle::writeFields(const Cloud<patchInteractSolidParticle>& c)
{
    particle::writeFields(c);

    const label np = c.size();
    const bool writeOnProc = c.size();

    IOField<scalar> d(c.fieldIOobject("d", IOobject::NO_READ), np);
    IOField<vector> U(c.fieldIOobject("U", IOobject::NO_READ), np);

    label i = 0;
    for (const auto& p : c)
    {
        d[i] = p.d_;
        U[i] = p.U_;
        ++i;
    }

    d.write(writeOnProc);
    U.write(writeOnProc);
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Foam::Ostream& os, const Foam::patchInteractSolidParticle& p)
{
    if (os.format() == IOstreamOption::ASCII)
    {
        os << static_cast<const particle&>(p) << token::SPACE << p.d_ << token::SPACE << p.U_;
    }
    else
    {
        os << static_cast<const particle&>(p);
        os.write(reinterpret_cast<const char*>(&p.d_), patchInteractSolidParticle::sizeofFields);
    }

    os.check(FUNCTION_NAME);
    return os;
}

// ************************************************************************* //
