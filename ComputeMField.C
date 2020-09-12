/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "ComputeMField.h"
#include "Node.h"
#include "SimParameters.h"
#include "HomePatch.h"
#include "Sequencer.h"
#include "HomePatch.h"
#include "ReductionMgr.h"
#include "CollectionMgr.h"
#include "BroadcastObject.h"
#include "Output.h"
#include "Controller.h"
#include "Broadcasts.h"
#include "Molecule.h"
#include "NamdOneTools.h"
#include "LdbCoordinator.h"
#include "Thread.h"
#include "Random.h"
#include "PatchMap.inl"
#include "ComputeMgr.h"
#include "ComputeGlobal.h"

ComputeMField::ComputeMField(ComputeID c, PatchID pid)
  : ComputeHomePatch(c,pid)
{

	reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);

}
/*			END OF FUNCTION ComputeMField		*/


ComputeMField::~ComputeMField()

{
	delete reduction;
}
/*			END OF FUNCTION ~ComputeMField		*/


void ComputeMField::doForce(FullAtom* p, Results* r) {

  SimParameters *simParams = Node::Object()->simParameters;
  Vector mField = simParams->mField;
  // Calculate the angular frequency in 1/fs.
  //TODO: Kosar, does frequency and phase  make sense in magnetic field? If not, change where approperate
  BigReal omega = TWOPI * simParams->mFieldFreq / 1000.;
  BigReal phi = PI/180.* simParams->mFieldPhase;
  BigReal t = patch->flags.step * simParams->dt; 
  Vector mField1 = cos(omega * t - phi) * mField;

  const int normalized = simParams->mFieldNormalized;
  if ( normalized ) {
    Lattice &l = homePatch->lattice;
    mField1 = Vector(l.a_r()*mField1, l.b_r()*mField1, l.c_r()*mField1);
  }

  Force *forces = r->f[Results::normal];
  BigReal energy = 0;
  Force extForce = 0.;
  Tensor extVirial;

  //  Loop through and check each atom
  for (int i=0; i<numAtoms; i++) {
    Force force = PDBVELFACTOR * (-p[i].charge) * cross( mField1, p[i].velocity);
    forces[i] += force;
    Position vpos = homePatch->lattice.reverse_transform(
		p[i].position, p[i].transform );
    energy -= force * (vpos - homePatch->lattice.origin()); //TODO: Kosar, energy preservation should hold

    printf(", velocity: [%f, %f, %f]\n\n",  p[i].velocity.x,  p[i].velocity.y,  p[i].velocity.z);
    printf(", VELOCITY: [%f, %f, %f]\n\n",  PDBVELFACTOR * p[i].velocity.x, PDBVELFACTOR * p[i].velocity.y,  PDBVELFACTOR * p[i].velocity.z);


    
    if ( ! normalized ) {
      extForce += force;
      extVirial += outer(force,vpos);
      }
}                                                                                                                                                            

    
  reduction->item(REDUCTION_MISC_ENERGY) += energy;
  if ( ! normalized ) {
    ADD_VECTOR_OBJECT(reduction,REDUCTION_EXT_FORCE_NORMAL,extForce);
    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NORMAL,extVirial);
  }
  reduction->submit();

}
/*			END OF FUNCTION force				*/
