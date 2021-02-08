{
    alphat = turbulence->nut()/Prt;
    alphat.correctBoundaryConditions();

    volScalarField alphaEff("alphaEff", turbulence->nu()/Pr + alphat);

    //get the grid spacings
    //const surfaceVectorField& deltas = mesh().delta();
    //Enthalpy update equation
    DH += cp*20*DiaCoeff/(0.0889*0.0666/(80*60)/runTime.deltaTValue())*(T - TRef);
    //DH += cp*1000*DiaCoeff/(0.0889*0.0666/(80*60)/runTime.deltaTValue())*(T - TRef);
    //DH += cp*10*(T - TRef);
    //Info<< "Coeff_1 = " << 0.2*DiaCoeff/(0.0889*0.0666/(80*60)/runTime.deltaTValue()) << endl;
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
    C_poro = C_mushy/7500*(1-f_l)*(1-f_l)/(f_l*f_l*f_l+eps_poro);

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(alphaEff, T)
     ==
        radiation->ST(rhoCpRef, T)
      + fvOptions(T)
      - fvc::ddt(DHcp)
     // + (-DHcp+DHcp_old)/runTime.deltaT()
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

forAll(DiaCoeff.internalField(),index)
    {
        DiaCoeff[index] = TEqn.diag()[index];
        //Info<< "DiaCoeff[" << index << "] = " << DiaCoeff[index] << endl;
    }
    

    TEqn.solve();
    //Info<< "dHdt_cp = " << (-DHcp+DHcp_old)/runTime.deltaT() << endl;

    radiation->correct();

    fvOptions.correct(T);

    rhok = 1.0 - beta*(T - TRef);

    //let DHcp_old = DHcp_new
    DHcp_old = DHcp;
}
