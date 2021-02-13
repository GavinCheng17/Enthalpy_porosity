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
     // - fvc::ddt(DHcp) //source term added to account for the latent energy
     // + (-DHcp+DHcp_old)/runTime.deltaT()
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    //store the diagonal coefficients and cell volumes for later use
    forAll(cell_vol.internalField(),index)
        {
            cell_vol[index] = mesh.V()[index];
            DiaCoeff[index] = TEqn.diag()[index];
        }

    double error;

    int iter = 0;

    do
    {
        Info<<  "Iteration = " << iter <<endl;
        //enthalpy update****************************
        
        /*//get the grid spacings
        const surfaceVectorField& deltas = mesh.delta();
        forAll(deltas, I)
        {
            Info<<  "dx[" << I << "] = " << deltas[I].x() <<endl;
            Info<<  "dy[" << I << "] = " << deltas[I].y() <<endl;
            Info<<  "dz[" << I << "] = " << deltas[I].z() <<endl;
        }
        */

        //Enthalpy update equation
        //DH += cp*20*DiaCoeff/(deltas.x()*deltas.y()/runTime.deltaTValue())*(T - TRef);
        DH = DHcp_old*cp + cp*omega_DH*DiaCoeff/(cell_vol/runTime.deltaTValue())*(T - TRef);
    
        //Ensure DH is bounded
        dimensionedScalar MaxDH = L_heat;
        dimensionedScalar MinDH("MinDH", dimEnergy/dimMass, 0.0);
        DH = max(min(DH, MaxDH), MinDH);
    
        //update the liquid fraction 
        f_l = DH/L_heat;
    
        //update DHcp
        DHcp = DH/cp;
        Info<<  "DH/cp = " << DH/cp <<endl;

        //update the porosity term
        C_poro = C_mushy/rho_ref*(1-2*f_l+f_l*f_l)/(f_l*f_l*f_l+eps_poro);
        //end of enthalpy update****************************
    
        Info<<  "Posi = 1"  <<endl;
        
        //reinitialise the temperature field
        T = T_old;

        solve
        (
            TEqn
         ==
        //    - fvc::ddt(DHcp) //source term added to account for the latent energy
             (-DHcp+DHcp_old)/runTime.deltaT()
        );

        Info<<  "Posi = 2"  <<endl;

        radiation->correct();

        fvOptions.correct(T);

        //rhok = 1.0 - beta*(T - TRef);

        //compute overall errors between the enthalpy-energy equation sweep 
        error = 0.0;
        forAll(DiaCoeff.internalField(),index)
        {
            error += fabs(T_star[index]-T[index]);
        }

        T_star = T;
        Info<<  "Enthalpy-energy sweep tolorance = " << error <<endl;
        iter++;
    }while (error >= 1e-7);

    //let DHcp_old = DHcp_new
    DHcp_old = DHcp;
    T_old = T;
    T_star = T;
}
