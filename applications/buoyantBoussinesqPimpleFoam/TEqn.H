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
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    //store diagonal coefficients of the T equation and the cell volumes to use for the next time iteration
forAll(DiaCoeff.internalField(),index)
    {
        DiaCoeff[index] = TEqn.diag()[index];
        cell_vol[index] = mesh.V()[index];
    }

    //enthalpy update****************************

    //Enthalpy update equation
    DH += cp*omega_DH*DiaCoeff/(cell_vol/runTime.deltaTValue())*(T - TRef);

    //Ensure DH is bounded
    dimensionedScalar MaxDH = L_heat;
    dimensionedScalar MinDH("MinDH", dimEnergy/dimMass, 0.0);
    DH = max(min(DH, MaxDH), MinDH);

    //update the liquid fraction 
    f_l = DH/L_heat;

    //update DHcp
    DHcp = DH/cp;

    //update the porosity term
    C_poro = C_mushy/rho_ref*(1-2*f_l+f_l*f_l)/(f_l*f_l*f_l+eps_poro);
    //end of enthalpy update****************************
    
    solve
    (
        TEqn
     ==
       - fvc::ddt(DHcp) //source term added to account for the latent energy
    );

    radiation->correct();

    fvOptions.correct(T);
}
