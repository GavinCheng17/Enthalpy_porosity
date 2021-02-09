singlePhaseTransportModel laminarTransport(U, phi);

// Thermal expansion coefficient [1/K]
dimensionedScalar beta
(
    "beta",
    dimless/dimTemperature,
    laminarTransport
);

// Reference temperature [K]
dimensionedScalar TRef("TRef", dimTemperature, laminarTransport);

// Laminar Prandtl number
dimensionedScalar Pr("Pr", dimless, laminarTransport);

// Turbulent Prandtl number
dimensionedScalar Prt("Prt", dimless, laminarTransport);

// Latent heat of fusion for tin
dimensionedScalar L_heat("L_heat", dimEnergy/dimMass, laminarTransport);

//Specific heat capacity
dimensionedScalar cp("cp", dimEnergy/dimMass/dimTemperature, laminarTransport);

//Mushy zone constant
dimensionedScalar C_mushy("C_mushy", dimDensity/dimTime, laminarTransport);

//Computational constant to avoid division by zero
dimensionedScalar eps_poro("eps_poro", dimless, laminarTransport);

//reference density kg/m^3
dimensionedScalar rho_ref("rho_ref", dimDensity, laminarTransport);

//dx*dy*dz used for coefficeint A_p^old in the enthalpy update equation
dimensionedScalar dxdydz("dxdydz", dimless, laminarTransport);

//relaxation factor in the enthalpy update equation
dimensionedScalar omega_DH("omega_DH", dimless, laminarTransport);
