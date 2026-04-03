#include "mapUtils.h"
#include "fvCFD.H"
#include "pisoControl.H"

int main(int argc, char* argv[])
{
    argList::noParallel();
    argList::addOption("text", "RenderText");
    argList::addOption("font", "FontPath");
    argList::addOption("fontsize", "FontSize");
    argList::addOption("scale", "ScaleFactor");

#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"

    const auto text = args.getOrDefault<word>("text", "Hello");
    const auto fontPath =
      args.getOrDefault<fileName>("font", "/usr/share/fonts/truetype/ubuntu/Ubuntu-M.ttf");
    auto fontSize = args.getOrDefault<float>("fontsize", 30.0f);

    auto scaleFactor = args.getOrDefault<scalar>("scale", 1.0);

    auto matrix = MapUtils::renderTextToMatrix(text, fontPath, fontSize);
    // flip the matrix upside down to match coordinate system
    std::reverse(matrix.begin(), matrix.end());
    Info << "Rendered text matrix (" << matrix.size() << "x" << matrix[0].size() << "):\n";

    // read the block sizes of the mesh
    auto meshDict = IOdictionary(IOobject(
      "blockMeshDict", runTime.system(), mesh, IOobject::MUST_READ_IF_MODIFIED, IOobject::NO_WRITE
    ));
    auto vec = MapUtils::parseCells(meshDict.lookup("blocks"));
    auto mappedMatrix = MapUtils::mapMatrix(std::move(matrix), {vec.x(), vec.y()});

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    const auto Nc = mesh.nCells();
    // Set the T value based on the rendered text matrix
    auto field = scalarField(Nc, 0.0);
    for (auto y = 0; y < static_cast<int>(mappedMatrix.size()); ++y)
    {
        for (auto x = 0; x < static_cast<int>(mappedMatrix[0].size()); ++x)
        {
            auto index = y * static_cast<int>(mappedMatrix[0].size()) + x;
            if (index < Nc)
            {
                field[index] =
                  static_cast<scalar>(mappedMatrix[y][x]) * mesh.V()[index] * scaleFactor;
            }
        }
    }

    // Create Fields
    auto transportProperties = IOdictionary(IOobject(
      "transportProperties",
      runTime.constant(),
      mesh,
      IOobject::MUST_READ_IF_MODIFIED,
      IOobject::NO_WRITE
    ));
    auto T = volScalarField(
      IOobject("T", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE), mesh
    );
    auto DT = dimensionedScalar("DT", dimLength * dimLength / dimTime, transportProperties);
    // Set the initial field values
    T.primitiveFieldRef() = field;
    T.correctBoundaryConditions();

    auto nu = dimensionedScalar("nu", dimViscosity, transportProperties);

    Info << "Reading field p\n" << endl;
    volScalarField p(
      IOobject("p", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh
    );

    Info << "Reading field U\n" << endl;
    volVectorField U(
      IOobject("U", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh
    );

    auto phi = surfaceScalarField(
      IOobject("phi", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
      fvc::flux(U)
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
    mesh.setFluxRequired(p.name());

    // Initilize
    auto cumulativeContErrIO = uniformDimensionedScalarField(
      IOobject(
        "cumulativeContErr",
        runTime.timeName(),
        "uniform",
        mesh.thisDb(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
      ),
      dimensionedScalar(word::null, dimless, Zero)
    );
    auto& cumulativeContErr = cumulativeContErrIO.value();

    auto piso = pisoControl(mesh);

    // Time loop
    while (runTime.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;

        // solving U
        {
            auto UEqn = (fvm::ddt(U) + fvm::div(phi, U) - fvm::laplacian(nu, U))();
            if (piso.momentumPredictor())
            {
                solve(UEqn == -fvc::grad(p));
            }

            // --- PISO loop
            while (piso.correct())
            {
                auto rAU = (1.0 / UEqn.A())();
                auto HbyA = (constrainHbyA(rAU * UEqn.H(), U, p))();
                auto phiHbyA = surfaceScalarField(
                  "phiHbyA", fvc::flux(HbyA) + fvc::interpolate(rAU) * fvc::ddtCorr(U, phi)
                );

                adjustPhi(phiHbyA, U, p);

                // Update the pressure BCs to ensure flux consistency
                constrainPressure(p, U, phiHbyA, rAU);

                // Non-orthogonal pressure corrector loop
                while (piso.correctNonOrthogonal())
                {
                    // Pressure corrector
                    auto pEqn = (fvm::laplacian(rAU, p) == fvc::div(phiHbyA))();
                    pEqn.setReference(pRefCell, pRefValue);
                    pEqn.solve(p.select(piso.finalInnerIter()));

                    if (piso.finalNonOrthogonalIter())
                    {
                        phi = phiHbyA - pEqn.flux();
                    }
                }

                {
                    auto contErr = volScalarField(fvc::div(phi));
                    auto sumLocalContErr =
                      runTime.deltaTValue() * mag(contErr)().weightedAverage(mesh.V()).value();
                    auto globalContErr =
                      runTime.deltaTValue() * contErr.weightedAverage(mesh.V()).value();
                    cumulativeContErr += globalContErr;

                    Info << "time step continuity errors : sum local = " << sumLocalContErr
                         << ", global = " << globalContErr << ", cumulative = " << cumulativeContErr
                         << endl;
                }

                U = HbyA - rAU * fvc::grad(p);
                U.correctBoundaryConditions();
            }
        }

        {
            // Construct the equation
            auto TEqn = (fvm::ddt(T) + fvm::div(phi, T) - fvm::laplacian(DT, T))();

            // Relax the equation
            TEqn.relax();

            // Solve the equation
            solve(TEqn);

            // Update boundary conditions
            T.correctBoundaryConditions();
        }

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    return 0;
}
