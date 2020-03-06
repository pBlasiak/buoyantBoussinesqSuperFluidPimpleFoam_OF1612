{
    //alphat = tnut()/Prt;
    alphat = turbulence->nut()/Prt;
    alphat.correctBoundaryConditions();

    //tmp<volScalarField> talphaEff
    //(
	//    superFluid->alpha()
    //);
	//alphaEff = talphaEff.ref();
	alphaEff = kappa/rho/cp;

	//k = superFluid->k();
	//const dimensionedScalar aaa("aaa", dimensionSet(3,1,-9,-2,0,0,0), 53760214065249.20);
	//const dimensionedScalar bbb("bbb", dimensionSet(3,1,-9,-1,0,0,0), 85119133182715.50);
	//volScalarField K = aaa*T - bbb;
	//Info<< "K = " << K << endl;
	//volScalarField magG = max(mag(G), dimensionedScalar("minDt", dimTemperature/dimLength, SMALL));
	//volScalarField onebyf = pow
	//	   (
	//	       max
	//	       (
	//	           K/magG/magG, 
	//	    	   dimensionedScalar("small", dimensionSet(3,3,-9,-3,0,0,0), SMALL)
	//	       ), 
	//	       1./3
	//	   );

	tmp<volScalarField> tGMpGrest
	(
		superFluid->GM()/rho/cp + fvc::div(kappa*superFluid->Grest(), "div(Grest)")/rho/cp
	);
	volScalarField& GMpGrest = tGMpGrest.ref();

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      //- fvm::laplacian(alphaEff, T)
      - 1./rho/cp*fvm::laplacian(kappa, T)
//      - fvm::Sp(fvc::div(phi),T) // Taken from Nima Sam (post nr 107)
     ==
        radiation->ST(rhoCpRef, T) // nie policzone rhoCpRef
	  + GMpGrest
//      + fvm::SuSp(GMpGrest/T, T)
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
