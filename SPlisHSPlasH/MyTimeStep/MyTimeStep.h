#ifndef __MyTimeStep_h__
#define __MyTimeStep_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "SPlisHSPlasH/SPHKernels.h"

namespace SPH
{
	class MyTimeStep : public TimeStep
	{
		private:
		unsigned int steps=0;
	public:
		MyTimeStep();
		virtual ~MyTimeStep(void);
		virtual void step() override;
		virtual void resize() override;
		virtual void reset() override;
	};
}

#endif
