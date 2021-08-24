{
    //alphat = tnut()/Prt;
    alphat = turbulence->nut()/Prt;
    alphat.correctBoundaryConditions();

    //tmp<volScalarField> talphaEff
    //(
	//    superFluid->alpha()
    //);
	//alphaEff = talphaEff.ref();
	alphaEff = kappa/rho/cp + alphat;
	kappa += alphat*rho*cp; // added

	tmp<volScalarField> tGMpGrest
	(
		superFluid->GM()/rho/cp + fvc::div(kappa*superFluid->Grest(), "div(Grest)")/rho/cp
	);
	volScalarField& GMpGrest = tGMpGrest.ref();
	
	//Info<< "GMpGrest = " << GMpGrest << endl;

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      //- fvm::laplacian(alphaEff, T)
      - 1./rho/cp*fvm::laplacian(kappa, T)
//      - fvm::Sp(fvc::div(phi),T) // Taken from Nima Sam (post nr 107)
     ==
        radiation->ST(rhoCpRef, T) // nie policzone rhoCpRef
//	  + GMpGrest
      - fvm::SuSp(-GMpGrest/T, T)
      + fvOptions(T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    radiation->correct();

    fvOptions.correct(T);
    T.correctBoundaryConditions();

	superFluid->correct();

	beta = superFluid->betaHe();
    rhok = 1.0 - beta*(T - TRef);
	cp = superFluid->cpHe();
	rho = superFluid->rhoHe();
	kappa = superFluid->alpha();
	Pr = superFluid->Pr(); // it is needed for alphat BC

    Info<< "min/max(T) = " << min(T).value() << ", "
        << max(T).value() <<endl;
}
