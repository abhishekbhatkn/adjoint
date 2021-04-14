{
    alphat = turbulence->nut()/Prt;
    alphat.correctBoundaryConditions();
    
    volScalarField alphaEff("alphaEff", alphaSolid+((turbulence->nu()/Pr)-alphaSolid)*eta*(1.0+penalty)/(eta+penalty)+alphat);

    fvScalarMatrix TEqn
    (
        fvm::div(phi, T)
      - fvm::laplacian(alphaEff, T)
     ==
        radiation->ST(rhoCpRef, T)
      + fvOptions(T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    radiation->correct();

    fvOptions.correct(T);

    rhok = 1.0 - beta*(T - TRef);
}
