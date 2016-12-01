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

void associate_ap2_clath_sigmaPBCCELL(int p1, int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, double bindrad, int *ncrosscom) {
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//NOTE: This function is was writing by Osman and Basilio. However, this is a modification of associate_freelegPBCELL.cpp fuction///
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	int prod = Rlist[mu][2];
	int iind = ihome[i1];
	int iind2 = ihome[i2];
	double Nv1;

	//get values of original associating complexes
	int c1 = bases[p1].mycomplex;
	int s1 = ind_com[c1].mysize;
	int c2 = bases[p2].mycomplex;
	int s2 = ind_com[c2].mysize;
	int newsize = s1 + s2;
	double trot;

	double *v1 = new double[3];
	double *v2 = new double[3];
	double *M = new double[9];
	double *Mneg = new double[9]; //same axis, negative angle
	double *u = new double[3];
	double *dtrans = new double[3];
	double *drev = new double[3];
	double cmpx=0.0, cmpy=0.0, cmpz=0.0;

	/*Now move protein*/
	double Dxtot = ind_com[c1].Dx + ind_com[c2].Dx;
	double Dytot = ind_com[c1].Dy + ind_com[c2].Dy;

	double Dztot = ind_com[c1].Dz + ind_com[c2].Dz;
	double tol = 1E-16;
	if (Dztot < tol)
		Dztot = 1; //otherwise you divide by zero

	//distance between the associating proteins interfaces!
	/*The sign needs to be switched to ensure the COM's get rotated
	 to the outside of interfaces final spot*/
	double dx = -bases[p2].x[iind2] + bases[p1].x[iind];
	double dy = -bases[p2].y[iind2] + bases[p1].y[iind];
	double dz = -bases[p2].z[iind2] + bases[p1].z[iind];

	dx -= plist.xboxl * round(dx / plist.xboxl);
	dy -= plist.yboxl * round(dy / plist.yboxl);
//	dz -= plist.zboxl * round(dz / plist.zboxl);

	double theta;

	///// THIS WAS WORKING CODE COULD BE USED
//	//distance to move to place interfaces, along vector v
	dtrans[0] = -dx * ind_com[c1].Dx / Dxtot;
	dtrans[1] = -dy * ind_com[c1].Dy / Dytot;
	dtrans[2] = -dz * ind_com[c1].Dz / Dztot;

	drev[0] = +dx * ind_com[c2].Dx / Dxtot;
	drev[1] = +dy * ind_com[c2].Dy / Dytot;
	drev[2] = +dz * ind_com[c2].Dz / Dztot;
//	double R2 = dx * dx + dy * dy + dz * dz;
//	double R1 = sqrt(R2);
//	//distance to move to place interfaces, along vector v
//	dtrans[0] = -dx * ind_com[c1].Dx / Dxtot* (bindrad / R1 - 1.0);
//	dtrans[1] = -dy * ind_com[c1].Dy / Dytot* (bindrad / R1 - 1.0);
//	dtrans[2] = -dz * ind_com[c1].Dz / Dztot* (bindrad / R1 - 1.0);
//
//	drev[0] = +dx * ind_com[c2].Dx / Dxtot* (bindrad / R1 - 1.0);
//	drev[1] = +dy * ind_com[c2].Dy / Dytot* (bindrad / R1 - 1.0);
//	drev[2] = +dz * ind_com[c2].Dz / Dztot* (bindrad / R1 - 1.0);
	if(bases[p1].protype == plist.pclath){//??? WHY JUST P1==? WHY NOT P2? It seems that if + ifelse (below) statement ensure that P1 or P2 are clath, if clath, change the center of mass????
		///////////////////////
		cout<<iind<<"EUREKA !!!"<<"BASILIO"<<endl;
		cout<<iind2<<"EUREKA !!!"<<"BASILIO"<<endl;
		/////////////////////////
		cmpx = bases[p1].x[iind-3]-bases[p1].xcom;
		cmpy = bases[p1].y[iind-3]-bases[p1].ycom;
		cmpz = bases[p1].z[iind-3]-bases[p1].zcom;
		cmpx -= plist.xboxl * round(cmpx / plist.xboxl);
		cmpy -= plist.yboxl * round(cmpy / plist.yboxl);
//BASILIO: Which values take iind? I Understand that iind is interfaces...so 1,2,3,4,5,6? what happen if iind take 1?
		v1[0] = bases[p1].xcom+cmpx/2 - bases[p1].x[iind];
		v1[1] = bases[p1].ycom+cmpy/2 - bases[p1].y[iind];
		v1[2] = bases[p1].zcom+cmpz/2 - bases[p1].z[iind];
		v1[0] -= plist.xboxl * round(v1[0] / plist.xboxl);
		v1[1] -= plist.yboxl * round(v1[1] / plist.yboxl);

		v2[0] = bases[p2].xcom - bases[p2].x[iind2];
		v2[1] = bases[p2].ycom - bases[p2].y[iind2];
		v2[2] = bases[p2].zcom - bases[p2].z[iind2];
		v2[0] -= plist.xboxl * round(v2[0] / plist.xboxl);
		v2[1] -= plist.yboxl * round(v2[1] / plist.yboxl);
	//	v2[2] -= plist.zboxl * round(v2[2] / plist.zboxl);
//		cout<<"new cm"<<"\n";
//		cout<<bases[p1].xcom+cmpx/2<<' '<<bases[p1].ycom+cmpy/2<<' '<<bases[p1].zcom+cmpz/2<<endl;
//		cout<<"old cm"<<"\n";
//		cout<<bases[p1].xcom<<' '<<bases[p1].ycom<<' '<<bases[p1].zcom<<endl;
//		cout<<"interface"<<"\n";
//		cout<<bases[p1].x[iind]<<' '<<bases[p1].y[iind]<<' '<<bases[p1].z[iind]<<endl;

	}else{

		////////////////////////////////
		cout<<iind<<"EUREKA !!!"<<"BASILIO"<<endl;
		cout<<iind2<<"EUREKA !!!"<<"BASILIO"<<endl;
		////////////////////////////////

		cmpx = bases[p2].x[iind2-3]-bases[p2].xcom;
		cmpy = bases[p2].y[iind2-3]-bases[p2].ycom;
		cmpz = bases[p2].z[iind2-3]-bases[p2].zcom;
		cmpx -= plist.xboxl * round(cmpx / plist.xboxl);
		cmpy -= plist.yboxl * round(cmpy / plist.yboxl);

		v2[0] = bases[p2].xcom+cmpx/2 - bases[p2].x[iind2];
		v2[1] = bases[p2].ycom+cmpy/2 - bases[p2].y[iind2];
		v2[2] = bases[p2].zcom+cmpz/2 - bases[p2].z[iind2];
		v2[0] -= plist.xboxl * round(v2[0] / plist.xboxl);
		v2[1] -= plist.yboxl * round(v2[1] / plist.yboxl);

		v1[0] = bases[p1].xcom - bases[p1].x[iind];
		v1[1] = bases[p1].ycom - bases[p1].y[iind];
		v1[2] = bases[p1].zcom - bases[p1].z[iind];
		v1[0] -= plist.xboxl * round(v1[0] / plist.xboxl);
		v1[1] -= plist.yboxl * round(v1[1] / plist.yboxl);
	//	v2[2] -= plist.zboxl * round(v2[2] / plist.zboxl);
//		cout<<"new cm"<<"\n";
//		cout<<bases[p2].xcom+cmpx/2<<' '<<bases[p2].ycom+cmpy/2<<' '<<bases[p2].zcom+cmpz/2<<endl;
//		cout<<"old cm"<<"\n";
//		cout<<bases[p2].xcom<<' '<<bases[p2].ycom<<' '<<bases[p2].zcom<<endl;
//		cout<<"interface"<<"\n";
//		cout<<bases[p2].x[iind2]<<' '<<bases[p2].y[iind2]<<' '<<bases[p2].z[iind2]<<endl;

	}

	/*Calculate rotation matrix*/

	//dotproduct(v2,v1, theta);
	double dp = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
	double l1 = v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2];
	double l2 = v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2];
	double dpnorm = dp / sqrt(l1 * l2);
	if (dpnorm > 1) {
		cout << "dp out of range! " << dpnorm << endl;
		theta = 0;
	} else if (dpnorm < -1) {
		cout << "dp out of range! " << dpnorm << endl;
		theta = M_PI;
	} else {
		theta = acos(dpnorm);
	}
	//cout<<theta<<"BASILIO"<<endl;
	//cout<<theta<<"EUREKA !!!"<<"BASILIO"<<endl;
	//now calculate the axis of rotation
	crossproduct(v1, v2, u); //u is the UNIT vector of the rotation axis
	//double trot = 0.5 * (M_PI - theta); //to push them 180 apart           ///////////BASILIO MODIFIED THIS
	/*for rotating v1 around u the negative trot */
//	ayudame(p1, p2, bases);
	////ADDED BY BASILIO///
	if(ind_com[c1].Dz!=0 && ind_com[c2].Dz!=0){	/*Update position of proteins in c2*/ // BOTH P1 AND P2 IN SOLUTION

		trot = 0.5*(M_PI - theta); //to push ONE 180. In this case clath

		//calc_Rmatrix(u, -trot,M);
		double ux = u[0];
		double uy = u[1];
		double uz = u[2];
		double cthet = cos(+trot);
		double sthet = sin(+trot);
		cout << "angle between: " << trot << " normal: " << ux << ' ' << uy << ' ' << uz << endl;
		M[0] = cthet + ux * ux * (1 - cthet);
		M[1] = ux * uy * (1 - cthet) - uz * sthet; //row 0, column 1, go across first!
		M[2] = ux * uz * (1 - cthet) + uy * sthet;
		M[3] = uy * ux * (1 - cthet) + uz * sthet;
		M[4] = cthet + uy * uy * (1 - cthet);
		M[5] = uy * uz * (1 - cthet) - ux * sthet;
		M[6] = uz * ux * (1 - cthet) - uy * sthet;
		M[7] = uz * uy * (1 - cthet) + ux * sthet;
		M[8] = cthet + uz * uz * (1 - cthet);

		Mneg[0] = M[0];
		Mneg[1] = ux * uy * (1 - cthet) + uz * sthet;
		Mneg[2] = ux * uz * (1 - cthet) - uy * sthet;
		Mneg[3] = uy * ux * (1 - cthet) - uz * sthet;
		Mneg[4] = M[4];
		Mneg[5] = uy * uz * (1 - cthet) + ux * sthet;
		Mneg[6] = uz * ux * (1 - cthet) + uy * sthet;
		Mneg[7] = uz * uy * (1 - cthet) - ux * sthet;
		Mneg[8] = M[8];
		rotate_and_translate_intPBCCELL(p1, c1, ind_com, bases, Mneg, dtrans, iind, plist);
		rotate_and_translate_intPBCCELL(p2, c2, ind_com, bases, M, drev, iind2, plist); //BASILIO: I AM MODIFYING THIS FUNCION

//		HERE WE ARE MISSING SIGMA DISPLACEMENT ALONG VECTOR V, TEST THIS
		//ayudame(p1, p2, bases);

		/*The sign ensures the COM's get rotated
		 to the outside of interfaces final spot*/

		if(bases[p1].protype == plist.pclath){

			cmpx = bases[p1].x[iind-3]-bases[p1].xcom;
			cmpy = bases[p1].y[iind-3]-bases[p1].ycom;
			cmpz = bases[p1].z[iind-3]-bases[p1].zcom;
			cmpx -= plist.xboxl * round(cmpx / plist.xboxl);
			cmpy -= plist.yboxl * round(cmpy / plist.yboxl);

			v1[0] = bases[p1].xcom+cmpx/2 - bases[p1].x[iind];
			v1[1] = bases[p1].ycom+cmpy/2 - bases[p1].y[iind];
			v1[2] = bases[p1].zcom+cmpz/2 - bases[p1].z[iind];
			v1[0] -= plist.xboxl * round(v1[0] / plist.xboxl);
			v1[1] -= plist.yboxl * round(v1[1] / plist.yboxl);

		}else{

			v1[0] = bases[p1].xcom - bases[p1].x[iind];
			v1[1] = bases[p1].ycom - bases[p1].y[iind];
			v1[2] = bases[p1].zcom - bases[p1].z[iind];
			v1[0] -= plist.xboxl * round(v1[0] / plist.xboxl);
			v1[1] -= plist.yboxl * round(v1[1] / plist.yboxl);
		//	v2[2] -= plist.zboxl * round(v2[2] / plist.zboxl);

		}
		Nv1 = sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);
		//distance to move to place interfaces, along vector v
		dtrans[0] = v1[0] * ind_com[c1].Dx / Dxtot * bindrad/Nv1;
		dtrans[1] = v1[1] * ind_com[c1].Dy / Dytot * bindrad/Nv1;
		dtrans[2] = v1[2] * ind_com[c1].Dz / Dztot * bindrad/Nv1;
		; //should be zero

		drev[0] = -v1[0] * ind_com[c2].Dx / Dxtot * bindrad/Nv1;
		drev[1] = -v1[1] * ind_com[c2].Dy / Dytot * bindrad/Nv1;
		drev[2] = -v1[2] * ind_com[c2].Dz / Dztot * bindrad/Nv1; //should be zero

		/*NOW WE ARE MOVING THE TWO PROTEINS TO A SEPARATION OF SIGMA*/
		//just translate if larger complexes are coming together
		translate_intPBCCELL(p1, c1, ind_com, bases, dtrans, plist);

		/*Same for complex2*********************************/
		//just translate
		translate_intPBCCELL(p2, c2, ind_com, bases, drev, plist);
//		ayudame(p1, p2, bases);

//////////////////////////////////////////////////////////////////////
		//THIS SECTION WAS WRITTEN BY BASILIO//
		//YOU CAN DELETE THIS SECTION BECAUSE IT IS NOT NECESARY, IT
		//WAS ADDED TO TEST THE ANGLE THETA///
		if(bases[p1].protype == plist.pclath){//??? WHY JUST P1==? WHY NOT P2? It seems that if + ifelse (below) statement ensure that P1 or P2 are clath, if clath, change the center of mass????

				cmpx = bases[p1].x[iind-3]-bases[p1].xcom;
				cmpy = bases[p1].y[iind-3]-bases[p1].ycom;
				cmpz = bases[p1].z[iind-3]-bases[p1].zcom;
				cmpx -= plist.xboxl * round(cmpx / plist.xboxl);
				cmpy -= plist.yboxl * round(cmpy / plist.yboxl);
		//BASILIO: Which values take iind? I Understand that iind is interfaces...so 1,2,3,4,5,6? what happen if iind take 1?
				v1[0] = bases[p1].xcom+cmpx/2 - bases[p1].x[iind];
				v1[1] = bases[p1].ycom+cmpy/2 - bases[p1].y[iind];
				v1[2] = bases[p1].zcom+cmpz/2 - bases[p1].z[iind];
				v1[0] -= plist.xboxl * round(v1[0] / plist.xboxl);
				v1[1] -= plist.yboxl * round(v1[1] / plist.yboxl);

				v2[0] = bases[p2].xcom - bases[p2].x[iind2];
				v2[1] = bases[p2].ycom - bases[p2].y[iind2];
				v2[2] = bases[p2].zcom - bases[p2].z[iind2];
				v2[0] -= plist.xboxl * round(v2[0] / plist.xboxl);
				v2[1] -= plist.yboxl * round(v2[1] / plist.yboxl);
			//	v2[2] -= plist.zboxl * round(v2[2] / plist.zboxl);
		//		cout<<"new cm"<<"\n";
		//		cout<<bases[p1].xcom+cmpx/2<<' '<<bases[p1].ycom+cmpy/2<<' '<<bases[p1].zcom+cmpz/2<<endl;
		//		cout<<"old cm"<<"\n";
		//		cout<<bases[p1].xcom<<' '<<bases[p1].ycom<<' '<<bases[p1].zcom<<endl;
		//		cout<<"interface"<<"\n";
		//		cout<<bases[p1].x[iind]<<' '<<bases[p1].y[iind]<<' '<<bases[p1].z[iind]<<endl;

			}else{

				cmpx = bases[p2].x[iind2-3]-bases[p2].xcom;
				cmpy = bases[p2].y[iind2-3]-bases[p2].ycom;
				cmpz = bases[p2].z[iind2-3]-bases[p2].zcom;
				cmpx -= plist.xboxl * round(cmpx / plist.xboxl);
				cmpy -= plist.yboxl * round(cmpy / plist.yboxl);

				v2[0] = bases[p2].xcom+cmpx/2 - bases[p2].x[iind2];
				v2[1] = bases[p2].ycom+cmpy/2 - bases[p2].y[iind2];
				v2[2] = bases[p2].zcom+cmpz/2 - bases[p2].z[iind2];
				v2[0] -= plist.xboxl * round(v2[0] / plist.xboxl);
				v2[1] -= plist.yboxl * round(v2[1] / plist.yboxl);

				v1[0] = bases[p1].xcom - bases[p1].x[iind];
				v1[1] = bases[p1].ycom - bases[p1].y[iind];
				v1[2] = bases[p1].zcom - bases[p1].z[iind];
				v1[0] -= plist.xboxl * round(v1[0] / plist.xboxl);
				v1[1] -= plist.yboxl * round(v1[1] / plist.yboxl);
			//	v2[2] -= plist.zboxl * round(v2[2] / plist.zboxl);
		//		cout<<"new cm"<<"\n";
		//		cout<<bases[p2].xcom+cmpx/2<<' '<<bases[p2].ycom+cmpy/2<<' '<<bases[p2].zcom+cmpz/2<<endl;
		//		cout<<"old cm"<<"\n";
		//		cout<<bases[p2].xcom<<' '<<bases[p2].ycom<<' '<<bases[p2].zcom<<endl;
		//		cout<<"interface"<<"\n";
		//		cout<<bases[p2].x[iind2]<<' '<<bases[p2].y[iind2]<<' '<<bases[p2].z[iind2]<<endl;

			}

			/*Calculate rotation matrix*/

			//dotproduct(v2,v1, theta);
			double dp = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
			double l1 = v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2];
			double l2 = v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2];
			double dpnorm = dp / sqrt(l1 * l2);
			if (dpnorm > 1) {
				cout << "dp out of range! " << dpnorm << endl;
				theta = 0;
			} else if (dpnorm < -1) {
				cout << "dp out of range! " << dpnorm << endl;
				theta = M_PI;
			} else {
				theta = acos(dpnorm);
			}
			//if (theta>0){
			if (theta<3.1 || theta > 3.16){
				cout<<theta<<"BASILIO"<<endl;
				cout<<theta<<"EUREKA !!!"<<"BASILIO"<<endl;
				ayudame(p1, p2, bases);
				cout<<"INCORRECT: My p1 is of type " << bases[p1].protype << endl;
			}else{
				ayudame(p1, p2, bases);
				cout<<"CORRECT: My p1 is of type " << bases[p1].protype << endl;
				//basilio
			}
		//////////////////////////////////////////////
		//DELETE TIL HERE WITHOUT PROBLEM////
		/////////////////////////////////////////

////////////////////////////////////////////////////////////////
		///COMENTED BY BASILIO 11-07-2016
	}else if(ind_com[c1].Dz==0 && ind_com[c2].Dz==0){ // everybody's on the membrane ///???? I BELIEVE IT HAVE TO BE bases[p1].Dz==0, NOT !=0
/// HERE FINISH COMMENT BY BASILIO 11-07-2016
		dx = -bases[p2].x[iind2] + bases[p1].x[iind];
		dy = -bases[p2].y[iind2] + bases[p1].y[iind];
		dz = -bases[p2].z[iind2] + bases[p1].z[iind];
		//BAISLIO: DELTAX
		dx -= plist.xboxl * round(dx / plist.xboxl);
		dy -= plist.yboxl * round(dy / plist.yboxl);
	//	dz -= plist.zboxl * round(dz / plist.zboxl);
		//boundary condition, inside box

		//distance to move to place interfaces, along vector v
		dtrans[0] = -dx * ind_com[c1].Dx / Dxtot;
		dtrans[1] = -dy * ind_com[c1].Dy / Dytot;
		dtrans[2] = 0.0;//-dz * ind_com[c1].Dz / Dztot;

		drev[0] = +dx * ind_com[c2].Dx / Dxtot;
		drev[1] = +dy * ind_com[c2].Dy / Dytot;
		drev[2] = 0.0;//+dz * ind_com[c2].Dz / Dztot;

		translate_intPBCCELL(p1, c1, ind_com, bases, dtrans, plist);
		translate_intPBCCELL(p2, c2, ind_com, bases, drev, plist);
// ???? I BELIEVE THAT NEXT SECTION SHOULD BE MODIFIED TO FIX PROBLEM CLATH + PIP2-AP2
	}else{//Clathrin in solution
		////what about clat-ap2 in solution interact with PI2?
		//ADDED BY BASILIO and osman
		double trot = (M_PI - theta); //0.5 * (M_PI - theta);  //to push them 180 apart
		/*for rotating v1 around u the negative trot */
		//???MODIFY THIS JUST FOR CLAT SOLVE THE PROBLEM. I mean the angle MPI-THETA JUST FOR CLATH

		//calc_Rmatrix(u, -trot,M);
		double ux = u[0];
		double uy = u[1];
		double uz = u[2];
		double cthet = cos(+trot);
		double sthet = sin(+trot);
		cout << "angle between: " << trot << " normal: " << ux << ' ' << uy << ' ' << uz << endl;
		M[0] = cthet + ux * ux * (1 - cthet);
		M[1] = ux * uy * (1 - cthet) - uz * sthet; //row 0, column 1, go across first!
		M[2] = ux * uz * (1 - cthet) + uy * sthet;
		M[3] = uy * ux * (1 - cthet) + uz * sthet;
		M[4] = cthet + uy * uy * (1 - cthet);
		M[5] = uy * uz * (1 - cthet) - ux * sthet;
		M[6] = uz * ux * (1 - cthet) - uy * sthet;
		M[7] = uz * uy * (1 - cthet) + ux * sthet;
		M[8] = cthet + uz * uz * (1 - cthet);
		//Mneg: Negative matrix??
		Mneg[0] = M[0];
		Mneg[1] = ux * uy * (1 - cthet) + uz * sthet;
		Mneg[2] = ux * uz * (1 - cthet) - uy * sthet;
		Mneg[3] = uy * ux * (1 - cthet) - uz * sthet;
		Mneg[4] = M[4];
		Mneg[5] = uy * uz * (1 - cthet) + ux * sthet;
		Mneg[6] = uz * ux * (1 - cthet) + uy * sthet;
		Mneg[7] = uz * uy * (1 - cthet) - ux * sthet;
		Mneg[8] = M[8];

		/*Update the position of proteins in complex one*/
		if(ind_com[c1].Dz!=0){//p1 is an clath molecule in the solution only rotate clathrin
			rotate_and_translate_intPBCCELL(p1, c1, ind_com, bases, Mneg, dtrans, iind, plist);
		}else{//if p2 is a clathrin molecule in the solution only rotate clathrin
		/*Update position of proteins in c2*/
			rotate_and_translate_intPBCCELL(p2, c2, ind_com, bases, M, drev, iind2, plist); //BASILIO: I AM MODIFYING THIS FUNCION
		}

	///////////////////////////////////////////////////////////////////////////
		dx = -bases[p2].x[iind2] + bases[p1].x[iind];
		dy = -bases[p2].y[iind2] + bases[p1].y[iind];
		dz = -bases[p2].z[iind2] + bases[p1].z[iind];
		//BAISLIO: DELTAX
		dx -= plist.xboxl * round(dx / plist.xboxl);
		dy -= plist.yboxl * round(dy / plist.yboxl);
	//	dz -= plist.zboxl * round(dz / plist.zboxl);
		//boundary condition, inside box

		//distance to move to place interfaces, along vector v
		dtrans[0] = -dx * ind_com[c1].Dx / Dxtot;
		dtrans[1] = -dy * ind_com[c1].Dy / Dytot;
		dtrans[2] = -dz * ind_com[c1].Dz / Dztot;

		drev[0] = +dx * ind_com[c2].Dx / Dxtot;
		drev[1] = +dy * ind_com[c2].Dy / Dytot;
		drev[2] = +dz * ind_com[c2].Dz / Dztot;

		translate_intPBCCELL(p1, c1, ind_com, bases, dtrans, plist);
		translate_intPBCCELL(p2, c2, ind_com, bases, drev, plist);
		dtrans[0]=0;
		dtrans[1]=0;
		dtrans[2]=bindrad;

		if(ind_com[c1].Dz!=0){///if p1 isnot Ap2(IS CLATHRIN) do this
			translate_intPBCCELL(p1, c1, ind_com, bases, dtrans, plist);
		}else{//if p2 isnot Ap2 do this
		/*Update position of proteins in c2*/
			translate_intPBCCELL(p2, c2, ind_com, bases, dtrans, plist);
		}
	///////////////////////////////////////////////////////////////////////////
//		ayudame(p1, p2, bases);
	}
	//// JUST CHECKING///////////////
	//////////////////////////////////////////////////////
	if(bases[p1].protype == plist.pclath){//??? WHY JUST P1==? WHY NOT P2? It seems that if + ifelse (below) statement ensure that P1 or P2 are clath, if clath, change the center of mass????

			cmpx = bases[p1].x[iind-3]-bases[p1].xcom;
			cmpy = bases[p1].y[iind-3]-bases[p1].ycom;
			cmpz = bases[p1].z[iind-3]-bases[p1].zcom;
			cmpx -= plist.xboxl * round(cmpx / plist.xboxl);
			cmpy -= plist.yboxl * round(cmpy / plist.yboxl);
	//BASILIO: Which values take iind? I Understand that iind is interfaces...so 1,2,3,4,5,6? what happen if iind take 1?
			v1[0] = bases[p1].xcom+cmpx/2 - bases[p1].x[iind];
			v1[1] = bases[p1].ycom+cmpy/2 - bases[p1].y[iind];
			v1[2] = bases[p1].zcom+cmpz/2 - bases[p1].z[iind];
			v1[0] -= plist.xboxl * round(v1[0] / plist.xboxl);
			v1[1] -= plist.yboxl * round(v1[1] / plist.yboxl);

			v2[0] = bases[p2].xcom - bases[p2].x[iind2];
			v2[1] = bases[p2].ycom - bases[p2].y[iind2];
			v2[2] = bases[p2].zcom - bases[p2].z[iind2];
			v2[0] -= plist.xboxl * round(v2[0] / plist.xboxl);
			v2[1] -= plist.yboxl * round(v2[1] / plist.yboxl);
		//	v2[2] -= plist.zboxl * round(v2[2] / plist.zboxl);
	//		cout<<"new cm"<<"\n";
	//		cout<<bases[p1].xcom+cmpx/2<<' '<<bases[p1].ycom+cmpy/2<<' '<<bases[p1].zcom+cmpz/2<<endl;
	//		cout<<"old cm"<<"\n";
	//		cout<<bases[p1].xcom<<' '<<bases[p1].ycom<<' '<<bases[p1].zcom<<endl;
	//		cout<<"interface"<<"\n";
	//		cout<<bases[p1].x[iind]<<' '<<bases[p1].y[iind]<<' '<<bases[p1].z[iind]<<endl;

		}else{

			cmpx = bases[p2].x[iind2-3]-bases[p2].xcom;
			cmpy = bases[p2].y[iind2-3]-bases[p2].ycom;
			cmpz = bases[p2].z[iind2-3]-bases[p2].zcom;
			cmpx -= plist.xboxl * round(cmpx / plist.xboxl);
			cmpy -= plist.yboxl * round(cmpy / plist.yboxl);

			v2[0] = bases[p2].xcom+cmpx/2 - bases[p2].x[iind2];
			v2[1] = bases[p2].ycom+cmpy/2 - bases[p2].y[iind2];
			v2[2] = bases[p2].zcom+cmpz/2 - bases[p2].z[iind2];
			v2[0] -= plist.xboxl * round(v2[0] / plist.xboxl);
			v2[1] -= plist.yboxl * round(v2[1] / plist.yboxl);

			v1[0] = bases[p1].xcom - bases[p1].x[iind];
			v1[1] = bases[p1].ycom - bases[p1].y[iind];
			v1[2] = bases[p1].zcom - bases[p1].z[iind];
			v1[0] -= plist.xboxl * round(v1[0] / plist.xboxl);
			v1[1] -= plist.yboxl * round(v1[1] / plist.yboxl);
		//	v2[2] -= plist.zboxl * round(v2[2] / plist.zboxl);
	//		cout<<"new cm"<<"\n";
	//		cout<<bases[p2].xcom+cmpx/2<<' '<<bases[p2].ycom+cmpy/2<<' '<<bases[p2].zcom+cmpz/2<<endl;
	//		cout<<"old cm"<<"\n";
	//		cout<<bases[p2].xcom<<' '<<bases[p2].ycom<<' '<<bases[p2].zcom<<endl;
	//		cout<<"interface"<<"\n";
	//		cout<<bases[p2].x[iind2]<<' '<<bases[p2].y[iind2]<<' '<<bases[p2].z[iind2]<<endl;

		}

		/*Calculate rotation matrix*/

		//dotproduct(v2,v1, theta);
		dp = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
		l1 = v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2];
		l2 = v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2];
		dpnorm = dp / sqrt(l1 * l2);
		if (dpnorm > 1) {
			cout << "dp out of range! " << dpnorm << endl;
			theta = 0;
		} else if (dpnorm < -1) {
			cout << "dp out of range! " << dpnorm << endl;
			theta = M_PI;
		} else {
			theta = acos(dpnorm);
		}
		if (theta<3.1 || theta > 3.16){
			cout<<theta<<"BASILIO 22 NOV 2016"<<endl;
			ayudame(p1, p2, bases);
			//basilio
		}
	/*Determine if you crashed two clathrins together
	 unless they are in the same complex, in which case you are closing a loop.
	 */
	int cancel;
	if (c1 != c2) {
		cancel = measure_overlap(c1, c2, ind_com, bases, plist.pclath);
	} else
		cancel = 0; //closing a loop, don't measure overlap becasue all proteins are in same complex

	if (cancel == 0) {
		/*Update the status, free and bound lists of these two proteins*/
		//change status of the interface
		bases[p1].istatus[iind] = prod;
		bases[p2].istatus[iind2] = prod;
		//no longer free
		int i;
		for (i = 0; i < bases[p1].nfree; i++) {
			if (bases[p1].freelist[i] == i1) {
				bases[p1].freelist[i] = bases[p1].freelist[bases[p1].nfree - 1]; //put last reaction in place of this one
				i = bases[p1].nfree;
			}
		}
		bases[p1].nfree -= 1;
		for (i = 0; i < bases[p2].nfree; i++) {
			if (bases[p2].freelist[i] == i2) {
				bases[p2].freelist[i] = bases[p2].freelist[bases[p2].nfree - 1]; //put last reaction in place of this one
				i = bases[p2].nfree;
			}
		}
		bases[p2].nfree -= 1;

		/*add this as a possible dissociation reaction to both proteins
		 so we know who they are bound to.
		 */
		bases[p1].bndlist[bases[p1].nbnd] = prod; //put this reaction as last possible for uni
		bases[p1].nbnd += 1; //now a dissociation reaction is possible

		bases[p2].bndlist[bases[p2].nbnd] = prod;
		bases[p2].nbnd += 1; //now a dissociation reaction is possible

		bases[p1].partner[iind] = p2;
		bases[p2].partner[iind2] = p1;

		/*add c2's proteins to c1's list, unless closing a loop*/
		/*Get new COM of complex*/
		/*Copy final complex into the spot of now gone complex c2 */
		int mp;
		int j;
		int tar;
		int flagbndry = 0;
		if (c1 == c2) {
			cout << "CLOSING A LOOP! " << endl;
			plist.nloop++;
			if (ind_com[c1].radR * 2.0 > plist.xboxl / 2.0)
				flagbndry = 1;
			update_one_com_onlyPBCCELL(c1, ind_com, bases, plist, flagbndry);

			update_radiusPBCCELL(c1, ind_com, bases, plist); //requires correct ind_com
			update_diffusion(c1, ind_com, bases);
////////////////////////////////////////////////////////////////////
//////////////////// MODIFY BY BASILIO 10/14/2016///////////////////
//			update_trans_diffusion(c1, ind_com, bases, plist.pretrans);
//////////////////////////////////////////////////////////////////////////
			//update_trans_diffusion(c1, ind_com, bases, plist.pretrans);
			update_rot_diffusion(c1, ind_com, bases, plist.prerot);
		} else {
			ind_com[c1].mysize = newsize;
			int t = s1;
			for (i = 0; i < s2; i++) {
				mp = ind_com[c2].plist[i];
				ind_com[c1].plist[t] = mp; //add these proteins to the first complex
				t++;
				bases[mp].mycomplex = c1;
			}
			if ((ind_com[c1].radR + ind_com[c2].radR) * 2.0 > plist.xboxl / 2.0)
				flagbndry = 1;
			update_one_com_onlyPBCCELL(c1, ind_com, bases, plist, flagbndry);
			update_radiusPBCCELL(c1, ind_com, bases, plist); //requires correct ind_com
			update_diffusion(c1, ind_com, bases);
////////////////////////////////////////////////////////////////////
//////////////////// MODIFY BY BASILIO 10/14/2016///////////////////
//			update_trans_diffusion(c1, ind_com, bases, plist.pretrans);
//////////////////////////////////////////////////////////////////////////
			update_rot_diffusion(c1, ind_com, bases, plist.prerot);

			plist.ntotalcomplex -= 1;
			tar = plist.ntotalcomplex;
			if (c2 != plist.ntotalcomplex) {

				/*otherwise, you are just deleting c2 entirely*/
				//copy element by element
				ind_com[c2].xcom = ind_com[tar].xcom;
				ind_com[c2].ycom = ind_com[tar].ycom;
				ind_com[c2].zcom = ind_com[tar].zcom;

				ind_com[c2].Dx = ind_com[tar].Dx;
				ind_com[c2].Dy = ind_com[tar].Dy;
				ind_com[c2].Dz = ind_com[tar].Dz;

				ind_com[c2].radR = ind_com[tar].radR;

				for (j = 0; j < ind_com[tar].mysize; j++) {
					ind_com[c2].plist[j] = ind_com[tar].plist[j];
				}
				ind_com[c2].mysize = ind_com[tar].mysize;

				for (j = 0; j < ind_com[c2].mysize; j++) {
					mp = ind_com[c2].plist[j];
					bases[mp].mycomplex = c2;
				}
			}
		}
	} else {
		/*IN THIS CASE< CLATHRINS CRASHED TOGETHER
		 un-bind them, allow to diffuse apart*/
		cout << "Unbind clathrins, crashed together ! " << endl;
		/*reverse previous rotation*/
		rotate_onlyPBCCELL(p1, c1, ind_com, bases, Mneg, plist);
		rotate_onlyPBCCELL(p2, c2, ind_com, bases, M, plist);

		dtrans[0] *= -1;
		dtrans[1] *= -1;
		dtrans[2] *= -1;

		drev[0] *= -1;
		drev[1] *= -1;
		drev[2] *= -1;

		/*Update the position of proteins in complex one*/
		rotate_and_translate_intPBCCELL(p1, c1, ind_com, bases, M, dtrans, iind, plist);

		/*Update position of proteins in c2*/
		rotate_and_translate_intPBCCELL(p2, c2, ind_com, bases, Mneg, drev, iind2, plist);

	}

	cout << "final complex com: " << ind_com[c1].xcom << ' ' << ind_com[c1].ycom << ' ' << ind_com[c1].zcom << " radius: " << ind_com[c1].radR << " Dr: " << ind_com[c1].Drx << " Dtrans: " << ind_com[c1].Dx << endl;

//	delete[] u1;
//	delete[] u2;
	delete[] v1;
	delete[] v2;
	delete[] u;
	delete[] M;
	delete[] Mneg;
	delete[] dtrans;
	delete[] drev;

}
