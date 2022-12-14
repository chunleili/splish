#ifndef __MyTimeStep_h__
#define __MyTimeStep_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "SimulationDataMyTimeStep.h"
#include "SPlisHSPlasH/SPHKernels.h"

namespace SPH
{
	class SimulationDataMyTimeStep;
	class MyTimeStep : public TimeStep
	{
	private:
		// SimulationDataMyTimeStep m_simulationData;
		unsigned int steps=0;
		void performNeighborhoodSearch();
		unsigned int m_counter;
	public:
		MyTimeStep();
		virtual ~MyTimeStep(void);
		virtual void step() override;
		virtual void resize() override;
		virtual void reset() override;
	};
}

#endif
