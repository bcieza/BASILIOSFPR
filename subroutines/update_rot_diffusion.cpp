#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include "reactions.h"
#include "vol_help.h"
#include "rand_gsl.h"
#include "Faddeeva.hh"
#include "utility_calls.h"
#include "vector_rot_calls.h"

void update_rot_diffusion(int c1, Complex *ind_com, Fullmol *bases, double prerot) {
// Previous version written by Prof. Johnson
//void update_rot_diffusion(int c1, Complex *ind_com, Fullmol *bases, double prerot) {
//	int i, p1, j;
//	int size = ind_com[c1].mysize;
//	double a = ind_com[c1].radR;
//	if (a == 0)
//		a = 1;
//	double a3 = a * a * a;
//	ind_com[c1].Drx = prerot / a3;
//	ind_com[c1].Dry = ind_com[c1].Drx;
//	ind_com[c1].Drz = ind_com[c1].Drx;
	/*Update the diffusion of a complex bases on sum of their radii,
	 defined via the Einstein-stokes equation, such that diffusions sum via their
	 inverses.
	 */
//////////////////////////////////////////////////////////////////////////////////
	int i, p1,j;
	int size = ind_com[c1].mysize;
	double Dxinv = 0;
	double Dyinv = 0;
	double Dzinv = 0;
	double inf = 1E500;
	p1 = ind_com[c1].plist[0];
	ind_com[c1].Drx = bases[p1].Drx;
	ind_com[c1].Dry = bases[p1].Dry;
	ind_com[c1].Drz = bases[p1].Drz;
	for (i = 0; i < size; i++) {

		p1 = ind_com[c1].plist[i];
		//Dxinv+=1.0/bases[p1].Dx;
		//Dyinv+=1.0/bases[p1].Dy;

		if (bases[p1].Drx != 0) {
			Dxinv += 1.0 / bases[p1].Drx;
		} else {
			Dxinv = inf;
			//      ind_com[c1].Dx=0;
		}
		if (bases[p1].Dry != 0) {
			Dyinv += 1.0 / bases[p1].Dry;
		} else {
			Dyinv = inf;
			//ind_com[c1].Dy=0;
		}
		if (bases[p1].Drz != 0) {
			Dzinv += 1.0 / bases[p1].Drz;
		} else {
			Dzinv = inf;
			//ind_com[c1].Dz=0;
		}
	}
	ind_com[c1].Drx = 1.0 / Dxinv;
	ind_com[c1].Dry = 1.0 / Dyinv;
	ind_com[c1].Drz = 1.0 / Dzinv;
	//explanation in 09 nov
}

