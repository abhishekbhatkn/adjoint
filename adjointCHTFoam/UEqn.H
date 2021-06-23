{
	forAll(eta,i){
		scalar valEta = eta[i];
		//scalar val = alpha_s+(alpha_f-alpha_s)*valEta*(1.0+penalty)/(valEta + penalty);
		porosity[i] = (porosity_s+(porosity_f-porosity_s)*valEta*(1.0+penalty)/(valEta + penalty));
	}

}
    // Solve the momentum equation

    MRF.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + MRF.DDt(U)
      + turbulence->divDevReff(U)
      + fvm::Sp(porosity, U)
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (simple.momentumPredictor())
    {
        solve
        (
            UEqn
          ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rhok)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
