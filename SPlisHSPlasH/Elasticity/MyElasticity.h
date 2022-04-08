#pragma once
#ifndef __MyElasticity_h__
#define __MyElasticity_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "ElasticityBase.h"

namespace SPH
{
	class MyElasticity : public ElasticityBase
	{
	protected:
		Real m_youngsModulus;
		Real m_poissonRatio;
		Vector3r m_fixedBoxMin;
		Vector3r m_fixedBoxMax;

		virtual void initParameters();
		void determineFixedParticles();

	public:
		MyElasticity(FluidModel* model);
		virtual ~MyElasticity(void);
		static int YOUNGS_MODULUS;
		static int POISSON_RATIO;
		static int FIXED_BOX_MIN;
		static int FIXED_BOX_MAX;

		static NonPressureForceBase* creator(FluidModel* model) { return new MyElasticity(model); }

		virtual void step();
		virtual void reset();
	};
}

#endif