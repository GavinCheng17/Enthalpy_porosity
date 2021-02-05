{
    alphat = turbulence->nut()/Prt;
    alphat.correctBoundaryConditions();


    //Enthalpy update equation
    DH += cp*runTime.deltaTValue()*(T - TRef);
    //Ensure DH is bounded
    dimensionedScalar MaxDH = L_heat;
    dimensionedScalar MinDH("MinDH", dimEnergy/dimMass, 0.0);
    DH = max(min(DH, MaxDH), MinDH);
    //update the liquid fraction 
    f_l = DH/L_heat;
    //update DHcp
    DHcp = DH/cp;
    //update the porosity term
    C_poro = -C_mushy*(1-f_l)*(1-f_l)/(f_l*f_l*f_l-eps_poro);

    volScalarField alphaEff("alphaEff", turbulence->nu()/Pr + alphat);

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(alphaEff, T)
     ==
        radiation->ST(rhoCpRef, T)
      + fvOptions(T)
      - (DHcp-DHcp_old)/runTime.deltaT()
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    radiation->correct();

    fvOptions.correct(T);

    rhok = 1.0 - beta*(T - TRef);

    //let DHcp_old = DHcp_new
    DHcp_old = DHcp;
}
