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

void rotate_and_translate_intPBCCELL(int p1, int c1, Complex *ind_com, Fullmol *bases, double *M, double *dtrans, int iind, Parms &plist) {

	/*First rotate all points in complex 1 around the com of protein p1, and translate them by vector d*/
	double *v = new double[3];
	double *v2 = new double[3];
	double *newpos = new double[3];
	int mp;
	int i, k;
	int s1 = ind_com[c1].mysize;
	//BASILIO: initial coordinates for p1
	double pivx = bases[p1].x[iind];
	double pivy = bases[p1].y[iind];
	double pivz = bases[p1].z[iind];

	for (i = 0; i < s1; i++) {
		mp = ind_com[c1].plist[i];
		for (k = 0; k < bases[mp].ninterface; k++) {
			// it begins to define coordinates for v vector?
			v[0] = bases[mp].x[k] - pivx;
			//v[0]=xcoordinates of interface 'k' minux pivx???
			v[1] = bases[mp].y[k] - pivy;
			v[2] = bases[mp].z[k] - pivz;
			v[0] -= plist.xboxl * round(v[0] / plist.xboxl);
			v[1] -= plist.yboxl * round(v[1] / plist.yboxl);
//			v[2] -= plist.zboxl * round(v[2] / plist.zboxl);

			rotate(v, M, v2); //includes the interface that will align  //BASILIO: I AM MODIFYING THIS FUNCION
			//BASILIO: I belive we have to comment rotate function() and change v2 by v in translate(v2,dtrans, newpos)
			translate(v2, dtrans, newpos);
			bases[mp].x[k] = pivx + newpos[0];
			bases[mp].y[k] = pivy + newpos[1];
			bases[mp].z[k] = pivz + newpos[2];
			bases[mp].x[k] -= plist.xboxl * round(bases[mp].x[k] / plist.xboxl);
			bases[mp].y[k] -= plist.yboxl * round(bases[mp].y[k] / plist.yboxl);
//			bases[mp].z[k] -= plist.zboxl * round(bases[mp].z[k] / plist.zboxl);

		}
		//rotate COM
		v[0] = bases[mp].xcom - pivx;
		v[1] = bases[mp].ycom - pivy;
		v[2] = bases[mp].zcom - pivz;
		v[0] -= plist.xboxl * round(v[0] / plist.xboxl);
		v[1] -= plist.yboxl * round(v[1] / plist.yboxl);
//		v[2] -= plist.zboxl * round(v[2] / plist.zboxl);

		rotate(v, M, v2); //includes the interface that will align
		//BASILIO: Again we can stop rotation to just translate the system
		translate(v2, dtrans, newpos);
		bases[mp].xcom = pivx + newpos[0];
		bases[mp].ycom = pivy + newpos[1];
		bases[mp].zcom = pivz + newpos[2];
		bases[mp].xcom -= plist.xboxl * round(bases[mp].xcom / plist.xboxl);
		bases[mp].ycom -= plist.yboxl * round(bases[mp].ycom / plist.yboxl);
//		bases[mp].zcom -= plist.zboxl * round(bases[mp].zcom / plist.zboxl);

	}
	delete[] v;
	delete[] v2;
	delete[] newpos;
}
