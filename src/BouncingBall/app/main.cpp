#include "fvCFD.H"
#include "patchInteractSolidParticleCloud.h"

int main(int argc, char* argv[])

{
#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
#include "readGravitationalAcceleration.H"

    // Create fields - these will be registered in the mesh database
    // and looked up by solidParticleCloud::move()
    Info << "\nReading field p\n" << endl;
    volScalarField p(
      IOobject("p", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh
    );

    Info << "Reading field U\n" << endl;
    volVectorField U(
      IOobject("U", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh
    );

    Info << "Reading field rho\n" << endl;
    volScalarField rho(
      IOobject("rho", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh
    );

    Info << "Reading field nu\n" << endl;
    volScalarField nu(
      IOobject("nu", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh
    );

    Info << "Creating particle cloud\n" << endl;
    patchInteractSolidParticleCloud cloud(mesh, "particleCloud");

    // Create particle with proper diameter and velocity
    // We'll create a temporary derived class to set the diameter

    auto particlepropertyDict = IOdictionary(IOobject(
      "particleProperties", runTime.constant(), mesh, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE
    ));

    auto particleDiameter = particlepropertyDict.get<List<scalar>>("diameter"); // 0.5 mm diameter
    auto particleVelocity = particlepropertyDict.get<List<vector>>("velocity"); // Initial velocity
    auto particlePosition = particlepropertyDict.get<List<vector>>("position"); // Initial position

    if (
      particleDiameter.size() != particlePosition.size() ||
      particleDiameter.size() != particleVelocity.size()
    )
    {
        FatalErrorInFunction << "Size mismatch in particle properties: "
                             << "diameter size = " << particleDiameter.size() << ", "
                             << "position size = " << particlePosition.size() << ", "
                             << "velocity size = " << particleVelocity.size() << exit(FatalError);
    }

    forAll(particleDiameter, i)
    {
        auto pos = particlePosition[i];
        auto vel = particleVelocity[i];
        auto d = particleDiameter[i];
        // Find the cell containing the particle position
        label celli = mesh.findCell(particlePosition[i]);
        if (celli < 0)
        {
            FatalErrorInFunction << "Particle position " << particlePosition
                                 << " is not inside the mesh domain." << exit(FatalError);
        }
        Info << "Adding particle at position " << particlePosition << " in cell " << celli
             << " with diameter " << particleDiameter << " m\n"
             << endl;

        // Create a temporary base particle with doLocate=true to properly initialize all
        // tetrahedral indices
        auto tempParticle = particle(mesh, pos, celli, -1, -1, true);

        // Now create the solid particle with the properly located coordinates
        cloud.addParticle(new patchInteractSolidParticle(
          mesh,
          tempParticle.coordinates(),
          tempParticle.cell(),
          tempParticle.tetFace(),
          tempParticle.tetPt(),
          d,
          vel
        ));
    }

    Info << cloud.size() << " particles created in the cloud.\n" << endl;

    while (runTime.run())
    {
        ++runTime;
        cloud.move(g);

        // Output particle information
        Info << "\nMoving particles in time = " << runTime.timeOutputValue() << " s" << nl << endl;
        forAllIter(patchInteractSolidParticleCloud, cloud, iter)
        {
            auto& particleI = iter();
            auto pos = particleI.position();
            auto pU = particleI.U();

            Info << "  position: " << pos << ", velocity: " << pU << nl;
        }
        runTime.write();
    }

    runTime.printExecutionTime(Info);
    Info << "\nEnd\n" << endl;

    return 0;
}
