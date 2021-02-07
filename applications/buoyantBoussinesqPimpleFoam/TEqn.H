{
    alphat = turbulence->nut()/Prt;
    alphat.correctBoundaryConditions();

    volScalarField alphaEff("alphaEff", turbulence->nu()/Pr + alphat);

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(alphaEff, T)
     ==
        radiation->ST(rhoCpRef, T)
      + fvOptions(T)
     // - fvc::ddt(DHcp)
      + (-DHcp+DHcp_old)/runTime.deltaT()
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    // create a scalar field to store diagonal coefficients 
    volScalarField DiaCoeff
(
    IOobject
    (
        "DiaCoeff",
        runTime.timeName(),
        mesh
    ), 
    f_l  //TEqn.diag()
);
forAll(DiaCoeff.internalField(),index)
    {
        DiaCoeff[index] = TEqn.diag()[index];
        //Info<< "DiaCoeff[" << index << "] = " << DiaCoeff[index] << endl;
    }
    //get the grid spacings
    //const surfaceVectorField& deltas = mesh().delta();
    //Enthalpy update equation
    DH += cp*0.2*DiaCoeff/(0.0889*0.0666/(80*60)/runTime.deltaTValue())*(T - TRef);
    //DH += cp*0.02*(T - TRef);
    //Info<< "Coeff_1 = " << 1/(0.0889*0.0666/(80*60)/runTime.deltaTValue())*(T - TRef) << endl;
    //Info<< "Coeff_2 = " << TEqn.diag() << endl;
    //Ensure DH is bounded
    dimensionedScalar MaxDH = L_heat;
    dimensionedScalar MinDH("MinDH", dimEnergy/dimMass, 0.0);
    DH = max(min(DH, MaxDH), MinDH);
    //Info<< "DH = " << DH << endl;
    //update the liquid fraction 
    f_l = DH/L_heat;
    //update DHcp
    DHcp = DH/cp;
    //update the porosity term
    C_poro = -C_mushy/7500*(1-f_l)*(1-f_l)/(f_l*f_l*f_l-eps_poro);

    TEqn.solve();
    //Info<< "Coeff_3 = " << TEqn.diag() << endl;

    radiation->correct();

    fvOptions.correct(T);

    rhok = 1.0 - beta*(T - TRef);

    //let DHcp_old = DHcp_new
    DHcp_old = DHcp;
}
