    // Solve the momentum equation

    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        fvm::ddt(U) + fvm::div(phi, U)
      + MRF.DDt(U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U)
      - fvm::Sp(C_poro,U) //porous source term
      - beta*g*(T-TRef) // nuoyant term
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                  - fvc::snGrad(p)
                )*mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
