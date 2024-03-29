/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fvc.H"
#include "KitamuraGradP.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace simplifiedSuperFluids
{
    defineTypeNameAndDebug(KitamuraGradP, 0);
    addToRunTimeSelectionTable(simplifiedSuperFluid, KitamuraGradP, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simplifiedSuperFluids::KitamuraGradP::KitamuraGradP
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    Kitamura(U, phi),
    KitamuraGradPCoeffs_(subDict(type() + "Coeffs")),
	p_(U.db().lookupObject<volScalarField>("p")),
	Grest_
    (
        IOobject
        (
            "Grest_",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
		U.mesh(),
		dimensionedVector("Grest", dimensionSet(0,-1,0,1,0,0,0), vector(0,0,0))
    ),
	GM_(KitamuraGradPCoeffs_.lookup("GM"))
{
//    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::simplifiedSuperFluids::KitamuraGradP::calcG()
{
	Grest_ = -fvc::grad(p_*rhoHe())/sHe()/rhoHe();
	G_ = fvc::grad(T_) + Grest_;
}

Foam::tmp<Foam::volScalarField> Foam::simplifiedSuperFluids::KitamuraGradP::GM() const
{
	if (GM_)
	{
		return 
		(
			   HeliumModelPtr_->AGM()*rhon_*rhos_*
		       pow
		       (
		           max
		           (
		               B_*magG_,
		         	    dimensionedScalar("small", dimensionSet(0,1,-1,0,0,0,0), SMALL)
		           ), 
		           4.0
		       )
		);
		//tmp<volScalarField>
		//(
		//    HeliumModelPtr_->AGM()*rhon_*rhos_*
		//    pow
		//    (
		//        max
		//        (
		//            B_*magG_,
		//      	    dimensionedScalar("small", dimensionSet(0,1,-1,0,0,0,0), SMALL)
		//        ), 
		//        4.0
		//    )
		//);
	} 
	else
	{
		return Foam::simplifiedSuperFluids::Kitamura::GM();
	}
}

Foam::tmp<Foam::volVectorField> Foam::simplifiedSuperFluids::KitamuraGradP::Grest() const
{
	return Grest_;
}

//Foam::tmp<Foam::volScalarField> Foam::simplifiedSuperFluids::KitamuraGradP::Grest() const
//{
//	return tmp<volScalarField>
//		   (
//		       fvc::div(Grest_, "div(Grest)") // it should be fvc::div(alphaEff*Grest_)
//		   );
//	return tmp<volScalarField>
//		   (
//		       new volScalarField
//			   (
//			       IOobject
//				   (
//				       T_.mesh().time().timeName(),
//					   T_.mesh(),
//					   IOobject::NO_READ,
//					   IOobject::NO_WRITE
//				   ),
//				   fvc::div(alpha()*Grest_, "div(Grest)") // it should be fvc::div(alphaEff*Grest_)
//			   )
//		   );
//}

bool Foam::simplifiedSuperFluids::KitamuraGradP::read()
{
    if (simplifiedSuperFluid::read())
    {
        KitamuraGradPCoeffs_ = subDict(type() + "Coeffs");

        KitamuraGradPCoeffs_.lookup("GM") >> GM_;

        return true;
    }
    else
    {
        return false;
    }
}
// ************************************************************************* //
