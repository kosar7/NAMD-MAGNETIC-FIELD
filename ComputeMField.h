/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTEMFIELD_H
#define COMPUTEMFIELD_H

#include "ComputeHomePatch.h"
#include "ReductionMgr.h"

class ComputeMField : public ComputeHomePatch
{

public:
	ComputeMField(ComputeID c, PatchID pid); 	//  Constructor
	virtual ~ComputeMField();			//  Destructor

	virtual void doForce(FullAtom* p, Results* r);

	SubmitReduction *reduction;

};

#endif







