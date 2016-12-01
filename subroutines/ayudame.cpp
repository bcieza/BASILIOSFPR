/*
 * ayudame.cpp
 *
 *  Created on: Oct 26, 2016
 *      Author: bcieza
 */

#include "reactions.h"
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

void ayudame(int p1, int p2, Fullmol *bases) {

	int nint = 0, n=0;
	cout<<"______"<<endl;
	cout<<bases[p1].xcom<<' '<<bases[p1].ycom<<' '<<bases[p1].zcom<<endl;
	if(p1<4101){
		//this is pip2
		nint = 2;
	}else{
		nint = 6;
	}
	for (n = 0; n < nint; n++) {
		cout << bases[p1].x[n] << ' ' << bases[p1].y[n] << ' ' << bases[p1].z[n] << endl;
	}

	cout<<bases[p2].xcom<<' '<<bases[p2].ycom<<' '<<bases[p2].zcom<<endl;
	if(p2<4101){
		//this is pip2
		nint = 2;
	}else{
		nint = 6;
	}
	for (n = 0; n < nint; n++) {
		cout << bases[p2].x[n] << ' ' << bases[p2].y[n] << ' ' << bases[p2].z[n] << endl;
	}
	cout<<"______"<<endl;


}



