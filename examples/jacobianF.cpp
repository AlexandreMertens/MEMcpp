#include <iostream>
#include <vector>

#include "TLorentzVector.h"
#include "TMath.h"

#include "jacobianF.h"

double computeJacobianF(const std::vector<TLorentzVector> p, const double sqrt_s){
	
	const double E1  = p.at(0).E();
	const double p1x = p.at(0).Px();
	const double p1y = p.at(0).Py();
	const double p1z = p.at(0).Pz();

	const double E2  = p.at(1).E();
	const double p2x = p.at(1).Px();
	const double p2y = p.at(1).Py();
	const double p2z = p.at(1).Pz();

	const double E3  = p.at(2).E();
	const double p3x = p.at(2).Px();
	const double p3y = p.at(2).Py();
	const double p3z = p.at(2).Pz();

	const double E4  = p.at(3).E();
	const double p4x = p.at(3).Px();
	const double p4y = p.at(3).Py();
	const double p4z = p.at(3).Pz();

    //double jac = (E4*(p1z*p2y*p3x - p1y*p2z*p3x - p1z*p2x*p3y + p1x*p2z*p3y + p1y*p2x*p3z - p1x*p2y*p3z) + E2*p1z*p3y*p4x - E1*p2z*p3y*p4x - E2*p1y*p3z*p4x + E1*p2y*p3z*p4x - E2*p1z*p3x*p4y + E1*p2z*p3x*p4y + E2*p1x*p3z*p4y - E1*p2x*p3z*p4y + (E2*p1y*p3x - E1*p2y*p3x - E2*p1x*p3y + E1*p2x*p3y)*p4z + E3*(-(p1z*p2y*p4x) + p1y*p2z*p4x + p1z*p2x*p4y - p1x*p2z*p4y - p1y*p2x*p4z + p1x*p2y*p4z));



    double inv_jac= (E4*(p1z*p2y*p3x - p1y*p2z*p3x - p1z*p2x*p3y + p1x*p2z*p3y + p1y*p2x*p3z - p1x*p2y*p3z) +  E2*p1z*p3y*p4x - E1*p2z*p3y*p4x - E2*p1y*p3z*p4x + E1*p2y*p3z*p4x - E2*p1z*p3x*p4y + E1*p2z*p3x*p4y +  E2*p1x*p3z*p4y - E1*p2x*p3z*p4y + (E2*p1y*p3x - E1*p2y*p3x - E2*p1x*p3y + E1*p2x*p3y)*p4z + E3*(-(p1z*p2y*p4x) + p1y*p2z*p4x + p1z*p2x*p4y - p1x*p2z*p4y - p1y*p2x*p4z + p1x*p2y*p4z));

    




        //if (inv_jac != jac) std::cout << "jac=" << jac  << " mw : " << inv_jac << std::endl;
	
	if(TMath::Abs(inv_jac) < JAC_MIN)
		return -1.;
	else
		return 1./( TMath::Abs(inv_jac) * 4.*16.*pow(TMath::Pi(),2.) );
}

