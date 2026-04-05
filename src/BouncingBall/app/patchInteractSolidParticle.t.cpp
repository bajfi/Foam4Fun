// Copyright (c) 2026 JackLee
//
// This software is released under the MIT License.
// https://opensource.org/licenses/MIT

#include "patchInteractSolidParticle.h"

namespace
{

using namespace Foam;

bool approxEqual(
  const scalar a, const scalar b, const scalar relTol = 1e-12, const scalar absTol = 1e-12
)
{
    const scalar diff = Foam::mag(a - b);
    if (diff <= absTol)
    {
        return true;
    }

    const scalar scale = Foam::max(Foam::mag(a), Foam::mag(b));
    return diff <= relTol * scale;
}

void expectVectorNear(
  const vector& actual,
  const vector& expected,
  std::string_view testName,
  std::size_t& passedChecks,
  const scalar relTol = 1e-12,
  const scalar absTol = 1e-12
)
{
    const bool ok = approxEqual(actual.x(), expected.x(), relTol, absTol) &&
                    approxEqual(actual.y(), expected.y(), relTol, absTol) &&
                    approxEqual(actual.z(), expected.z(), relTol, absTol);

    if (!ok)
    {
        const auto vecStr = [](const vector& v)
        {
            return "(" + std::to_string(v.x()) + ", " + std::to_string(v.y()) + ", " +
                   std::to_string(v.z()) + ")";
        };

        Info << "[FAIL] " << testName << '\n'
             << "  expected: " << vecStr(expected) << '\n'
             << "  actual  : " << vecStr(actual) << nl;
        std::exit(1);
    }

    ++passedChecks;
}

vector manualDragUpdate(
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

    scalar reFunc = 1.0;
    const scalar Re = magUr * d / nuc;
    if (Re > 0.01)
    {
        reFunc += 0.15 * std::pow(Re, 0.687);
    }

    const scalar Dc = (24.0 * nuc / d) * reFunc * (3.0 / 4.0) * (rhoc / (d * rhop));
    return (Up + dt * (Dc * Uc + (1.0 - rhoc / rhop) * g)) / (1.0 + dt * Dc);
}

} // namespace

int main()
{
    using namespace Foam;

    std::size_t passedChecks = 0;

    const patchInteractPhysicalProperties props{1000.0, 0.6, 0.1};

    // 1) dt=0 should preserve velocity exactly.
    {
        const vector Up(0.4, -0.3, 0.2);
        const vector Uc(1.0, 2.0, -1.0);
        const vector g(0.0, -9.81, 0.0);

        const vector actual = patchInteractSolidParticle::computeVelocityAfterDrag(
          Up, 1e-3, 0.0, 1.2, Uc, 1e-6, g, props
        );

        expectVectorNear(actual, Up, "computeVelocityAfterDrag_dt_zero", passedChecks);
    }

    // 2) Low-Re branch (Re <= 0.01).
    {
        const vector Up(0.0, 0.0, 0.0);
        const vector Uc(0.0, 0.0, 0.0); // magUr = 0 => Re = 0
        const vector g(0.0, -9.81, 0.0);

        const vector expected = manualDragUpdate(Up, 1e-3, 0.1, 1.0, Uc, 1e-6, g, props);
        const vector actual = patchInteractSolidParticle::computeVelocityAfterDrag(
          Up, 1e-3, 0.1, 1.0, Uc, 1e-6, g, props
        );

        expectVectorNear(actual, expected, "computeVelocityAfterDrag_low_Re", passedChecks);
    }

    // 3) High-Re branch (Re > 0.01).
    {
        const vector Up(0.0, 0.0, 0.0);
        const vector Uc(1.0, 0.0, 0.0); // high slip velocity
        const vector g(0.0, -9.81, 0.0);

        const vector expected = manualDragUpdate(Up, 1e-3, 0.05, 1.2, Uc, 1e-6, g, props);
        const vector actual = patchInteractSolidParticle::computeVelocityAfterDrag(
          Up, 1e-3, 0.05, 1.2, Uc, 1e-6, g, props
        );

        expectVectorNear(actual, expected, "computeVelocityAfterDrag_high_Re", passedChecks);
    }

    // 4) Wall interaction when moving into wall normal (Un > 0): bounce + tangential damping.
    {
        const vector Up(1.0, 2.0, 0.0);
        const vector wallNormal(0.0, 1.0, 0.0);

        const vector actual =
          patchInteractSolidParticle::computeVelocityAfterWallInteraction(Up, wallNormal, props);

        const vector expected(1.0 - props.mu, -2.0 * props.e, 0.0);
        expectVectorNear(
          actual, expected, "computeVelocityAfterWallInteraction_rebound", passedChecks
        );
    }

    // 5) Wall interaction when moving away from wall normal (Un <= 0): no bounce, tangential
    // damping only.
    {
        const vector Up(1.0, -2.0, 0.0);
        const vector wallNormal(0.0, 1.0, 0.0);

        const vector actual =
          patchInteractSolidParticle::computeVelocityAfterWallInteraction(Up, wallNormal, props);

        const vector expected(1.0 - props.mu, -2.0, 0.0);
        expectVectorNear(
          actual, expected, "computeVelocityAfterWallInteraction_no_rebound", passedChecks
        );
    }

    Info << "[PASS] patchInteractSolidParticle tests passed: " << passedChecks << " checks" << nl;
    return 0;
}
