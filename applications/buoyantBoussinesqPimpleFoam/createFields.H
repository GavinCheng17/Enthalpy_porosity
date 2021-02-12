Info<< "Reading thermophysical properties\n" << endl;

Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Read the nodal value of latent heat (Delat H)
Info<< "Reading field DH\n" << endl;
volScalarField DH
(
    IOobject
    (
        "DH",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
#include "createPhi.H"

#include "readTransportProperties.H"

Info<< "Creating turbulence model\n" << endl;
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

// Kinematic density for buoyancy force
volScalarField rhok
(
    IOobject
    (
        "rhok",
        runTime.timeName(),
        mesh
    ),
    1.0 - beta*(T - TRef)
);

// kinematic turbulent thermal thermal conductivity m2/s
Info<< "Reading field alphat\n" << endl;
volScalarField alphat
(
    IOobject
    (
        "alphat",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// liquid fraction
volScalarField f_l
(
    IOobject
    (
        "f_l",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    DH/L_heat
);

// porous term coefficient C_poro
volScalarField C_poro
(
    IOobject
    (
        "C_poro",
        runTime.timeName(),
        mesh
    ),
    C_mushy/rho_ref*(1-2*f_l+f_l*f_l)/(f_l*f_l*f_l+eps_poro)
);

/*
// create a scalar field to store the values of DH/cp at previous time step
volScalarField DHcp_old
(
    IOobject
    (
        "DHcp_old",
        runTime.timeName(),
        mesh
    ),
    DH/cp
);
*/

// create a scalar field to store the values of DH/cp
volScalarField DHcp
(
    IOobject
    (
        "DHcp",
        runTime.timeName(),
        mesh
    ),
    DH/cp
);

// create a scalar field to store diagonal coefficients at current time step Ap
    volScalarField DiaCoeff
(
    IOobject
    (
        "DiaCoeff",
        runTime.timeName(),
        mesh
    ), 
    f_l*0.0  //initial value as a non-dimensional number
);

// create a scalar field to store cell volumes
    volScalarField cell_vol
(
    IOobject
    (
        "cell_vol",
        runTime.timeName(),
        mesh
    ), 
    f_l*1.0  //initial value as a non-dimensional number
);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rhok*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
}

mesh.setFluxRequired(p_rgh.name());

#include "createMRF.H"
#include "createIncompressibleRadiationModel.H"
#include "createFvOptions.H"
