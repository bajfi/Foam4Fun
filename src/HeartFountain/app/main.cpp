// Copyright (c) 2026 JackLee
//
// This software is released under the MIT License.
// https://opensource.org/licenses/MIT

#include "fvCFD.H"
#include "pimpleControl.H"
#include "basicKinematicCollidingCloud.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char* argv[])
{
#include "addCheckCaseOptions.H"
#include "setRootCaseLists.H"
#include "createTime.H"
#include "createMesh.H"
#include "createTimeControls.H"
#include "readGravitationalAcceleration.H"

    pimpleControl pimple(mesh);

    Info << "Reading transportProperties\n" << endl;

    IOdictionary transportProperties(IOobject(
      "transportProperties",
      runTime.constant(),
      mesh,
      IOobject::MUST_READ_IF_MODIFIED,
      IOobject::NO_WRITE
    ));

    Info << "Reading field p\n" << endl;
    volScalarField p(
      IOobject("p", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, pimple.dict(), pRefCell, pRefValue);
    mesh.setFluxRequired(p.name());

    Info << "Reading field U\n" << endl;
    volVectorField U(
      IOobject("U", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh
    );

    Info << "Reading/calculating continuous-phase face flux field phi\n" << endl;

    surfaceScalarField phi(
      IOobject("phi", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
      linearInterpolate(U) & mesh.Sf()
    );

    Info << "Creating turbulence model\n" << endl;
    singlePhaseTransportModel continuousPhaseTransport(U, phi);

    dimensionedScalar rhocValue("rho", dimDensity, continuousPhaseTransport);
    volScalarField rhoc(
      IOobject(rhocValue.name(), runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
      mesh,
      rhocValue
    );

    volScalarField muc(
      IOobject("mu", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
      rhoc * continuousPhaseTransport.nu()
    );

    Info << "Creating field alphac\n" << endl;
    // alphac must be constructed before the cloud
    // so that the drag-models can find it
    volScalarField alphac(
      IOobject("alpha", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
      mesh,
      dimensionedScalar(dimless, Zero)
    );

    Info << "Constructing kinematicCloud " << endl;
    basicKinematicCollidingCloud kinematicCloud("kinematicCloud", rhoc, U, muc, g);

    // Particle fraction upper limit
    scalar alphacMin(
      1.0 -
      (kinematicCloud.particleProperties().subDict("constantProperties").get<scalar>("alphaMax"))
    );

    // Update alphac from the particle locations
    alphac = max(1.0 - kinematicCloud.theta(), alphacMin);
    alphac.correctBoundaryConditions();

    surfaceScalarField alphacf("alphacf", fvc::interpolate(alphac));
    surfaceScalarField alphaPhic("alphaPhi", alphacf * phi);

    auto continuousPhaseTurbulence(
      incompressible::turbulenceModel::New(U, phi, continuousPhaseTransport)
    );

    auto freeze_flow = transportProperties.lookupOrDefault<bool>("freezeFlow", true);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting time loop\n" << endl;

    scalar CoNum{};

    while (runTime.run())
    {
#include "readTimeControls.H"
#include "createPhi.H"
#include "setDeltaT.H"

        runTime++;
        Info << "Time = " << runTime.timeName() << nl << nl;

        continuousPhaseTransport.correct();
        muc = rhoc * continuousPhaseTransport.nu();

        Info << "Evolving cloud\n" << kinematicCloud.name() << nl;
        kinematicCloud.evolve();

        // Update continuous phase volume fraction field
        alphac = max(1.0 - kinematicCloud.theta(), alphacMin);
        alphac.correctBoundaryConditions();
        alphacf = fvc::interpolate(alphac);
        alphaPhic = alphacf * phi;

        Info << "Continuous phase Average volume fraction = "
             << alphac.weightedAverage(mesh.Vsc()).value()
             << "  Min(alphac) = " << min(alphac).value()
             << "  Max(alphac) = " << max(alphac).value() << nl;

        fvVectorMatrix cloudSU(kinematicCloud.SU(U));
        volVectorField cloudVolSUSu(
          IOobject("cloudVolSUSu", runTime.timeName(), mesh),
          mesh,
          dimensionedVector(cloudSU.dimensions() / dimVolume, Zero),
          zeroGradientFvPatchVectorField::typeName
        );

        cloudVolSUSu.primitiveFieldRef() = -cloudSU.source() / mesh.V();
        cloudVolSUSu.correctBoundaryConditions();
        cloudSU.source() = vector::zero;

        while (!freeze_flow && pimple.loop())
        {
            // solve the momentum equation for the continuous phase
            fvVectorMatrix UcEqn(
              fvm::ddt(alphac, U) + fvm::div(alphaPhic, U) -
                fvm::Sp(fvc::ddt(alphac) + fvc::div(alphaPhic), U) +
                continuousPhaseTurbulence->divDevRhoReff(U) ==
              (1.0 / rhoc) * cloudSU
            );

            UcEqn.relax();

            volScalarField rAUc(1.0 / UcEqn.A());
            surfaceScalarField rAUcf("Dp", fvc::interpolate(rAUc));

            surfaceScalarField phiForces(
              fvc::flux(rAUc * cloudVolSUSu / rhoc) + rAUcf * (g & mesh.Sf())
            );

            if (pimple.momentumPredictor())
            {
                solve(UcEqn == fvc::reconstruct(phiForces / rAUcf - fvc::snGrad(p) * mesh.magSf()));
            }

            while (pimple.correct())
            {
                volVectorField HbyA(constrainHbyA(rAUc * UcEqn.H(), U, p));

                surfaceScalarField phiHbyA(
                  "phiHbyA", (fvc::flux(HbyA) + alphacf * rAUcf * fvc::ddtCorr(U, phi))
                );

                if (p.needReference())
                {
                    adjustPhi(phiHbyA, U, p);
                }

                phiHbyA += phiForces;

                // Update the pressure BCs to ensure flux consistency
                constrainPressure(p, U, phiHbyA, rAUcf);

                // Non-orthogonal pressure corrector loop
                while (pimple.correctNonOrthogonal())
                {
                    fvScalarMatrix pEqn(
                      fvm::laplacian(alphacf * rAUcf, p) ==
                      fvc::ddt(alphac) + fvc::div(alphacf * phiHbyA)
                    );

                    pEqn.setReference(pRefCell, pRefValue);

                    pEqn.solve(p.select(pimple.finalInnerIter()));

                    if (pimple.finalNonOrthogonalIter())
                    {
                        phi = phiHbyA - pEqn.flux() / alphacf;

                        p.relax();

                        U = HbyA +
                            rAUc * fvc::reconstruct((phiForces - pEqn.flux() / alphacf) / rAUcf);
                        U.correctBoundaryConditions();
                    }
                }
            }

            if (pimple.turbCorr())
            {
                continuousPhaseTurbulence->correct();
            }
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
