#include "reactions.h"



void perform_association_sigma_com(int p1, Fullmol *bases, Complex *ind_com, Parms &plist, int **crosspart, int **cross_rxn, int **crossint, int **Rlist, int *i_home, double **probvec, int *Ncoup, int **mycoupled, int *p_home, int ci1, int ci2, int it, int *movestat, int *ncross, int bindrad, int *ncrosscom, double **traj)
{
  int p2 = crosspart[p1][ci1];//?????? WHY P2 AS FUNCTION OF P1
  int rxn1 = cross_rxn[p1][ci1];
  
  /*Associate proteins, move them to contact and update their free and bound lists//BASILIO: What are the free and bound list?
    For binding to lipids (they have D=0), move to contact but align vertically.
   */
  
  /*test for being in same complex to avoid moving binding interfaces too far*/
  //flag2=same_complex_test(p1, p2, bases, crossint, i_home, bindrad, rxn1, ci1, ci2);
  cout << "Associate between proteins: p1 " << p1 << ' ' << p2 << " interfaces: " << crossint[p1][ci1] << ' ' << crossint[p2][ci2] << " reaction: " << rxn1 << " pact; " << probvec[p1][ci1] <<  " iter: " << it << endl;
//BASILIO: crossint[p1][ci1]  interface ci1 of protein p1...
  /*for self binding, use freeleg for clathrin type proteins.*/
  //BASILIO REVIEW IT (OSAMN SAID DO IT!!!)

  ////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////// MODIFY BY BASILIO ////////////////////////////////////
  if (bases[p1].protype == plist.pclath && bases[p2].protype == plist.pclath){
	  if(ind_com[bases[p1].mycomplex].Dz!=0 && ind_com[bases[p2].mycomplex].Dz!=0) {
		  associate_freelegPBCCELL(p1, p2, rxn1, crossint[p1][ci1], crossint[p2][ci2], bases, Rlist, i_home, ind_com, plist);// 3d versus 3d clat clat
	  } else if (ind_com[bases[p1].mycomplex].Dz==0 && ind_com[bases[p2].mycomplex].Dz==0) {
		  associate_freelegPBCCELL(p1, p2, rxn1, crossint[p1][ci1], crossint[p2][ci2], bases, Rlist, i_home, ind_com, plist);// 2d versus 2d clat clat
	  } else {
		  associate_2dclat_3dclat_lipid_sigmaPBCCELL(p1, p2, rxn1, crossint[p1][ci1], crossint[p2][ci2], bases, Rlist, i_home, ind_com, plist, bindrad, ncrosscom); // 3d vs 2d and 2d vs 3d clat clat
	  }
  }else if (bases[p1].protype == plist.pclath || bases[p2].protype == plist.pclath) { //special ap2 clath
	  associate_ap2_clath_sigmaPBCCELL(p1, p2, rxn1, crossint[p1][ci1], crossint[p2][ci2], bases, Rlist, i_home, ind_com, plist, bindrad, ncrosscom);
  }else if(bases[p1].Dz==0 || bases[p2].Dz==0){//at least one partner is in 2D always.//??? This is just for PI2 ? Clath can be in the membrane too...
	  //???What Dz==0 means? No diffusion in z axes right? So this work for PI2, PI2-AP2 complex and P12-AP2-CLATH complex?
	  if(bases[p1].Dz==0 && bases[p2].Dz==0){
		  associate_zsigmaPBCCELL(p1, p2, rxn1, crossint[p1][ci1], crossint[p2][ci2], bases, Rlist, i_home, ind_com, plist, bindrad, ncrosscom);//2D association
	  }else{
		  associate_lipid_sigmaPBCCELLAP2(p1, p2, rxn1, crossint[p1][ci1], crossint[p2][ci2], bases, Rlist, i_home, ind_com, plist, bindrad, ncrosscom);//this comes from associate_freelegPBCCELL
	  }
	  // do you need the next section?
  }else
	  associate_zsigmaPBCCELL(p1, p2, rxn1, crossint[p1][ci1], crossint[p2][ci2], bases, Rlist, i_home, ind_com, plist, bindrad, ncrosscom); // when I modified by associate_ap2_clath_sigmaPBCCEL ...
  // ap2 interact not bertically with pip2
  ///////////////////// MODIFICATION END HERE /////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////
    
  ap2_coupledrxn_add(rxn1, Ncoup, mycoupled, bases, Rlist, i_home, p_home, p1, p2);
  reflect_complex_rad_rotCELL(p1, bases, ind_com, plist.xboxl, plist.yboxl, plist.zboxl);
  
  /*Remove p1 and p2 from the list of potential reaction partners so they don't try to
    associate again in this turn. They will also not try to avoid overlap anymore, and proteins
    close by will not try to avoid overlapping them this turn. 
  */
  //  remove_reaction_all_ncom(p1, p2, ncross, crosspart, probvec,  cross_rxn, crossint, ncrosscom, bases);
  /*
    To allow proteins that just associated to continue avoiding overlap, set probabilities to zero, and trajectories to zero.
  */
  traj[bases[p1].mycomplex][0]=0;
  traj[bases[p1].mycomplex][1]=0;
  traj[bases[p1].mycomplex][2]=0;
  /*Set probability to zero */
  remove_one_prob_all(p1, ncross, crosspart, probvec, cross_rxn, crossint);
  /*Set probability to zero */
  remove_one_prob_all(p2, ncross, crosspart, probvec, cross_rxn, crossint);
  ncross[p1]=-1;
  ncross[p2]=-1;
  ncrosscom[bases[p1].mycomplex]=-1;//you won't avoid them, but they should avoid you. you won't move.
  ncrosscom[bases[p2].mycomplex]=-1;
  
  /*Since these proteins have moved to associate and taken their complex with them,
    don't allow any proteins in their complex to move again.
  */
  set_movestat_zero(p1, bases, ind_com, movestat);
  
}
