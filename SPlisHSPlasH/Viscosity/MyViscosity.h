#pragma once
#ifndef __MyViscosity_h__
#define __MyViscosity_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "ViscosityBase.h"

namespace SPH
{
	class MyViscosity : public ViscosityBase
	{
	protected:
		virtual void initParameters();

		std::vector<Vector3r> m_vDiff;

		/** This function is called after the simulation scene is loaded and all
		* parameters are initialized. While reading a scene file several parameters
		* can change. The deferred init function should initialize all values which
		* depend on these parameters.
		*/
		virtual void deferredInit();

		void initValues();


	public:
		MyViscosity(FluidModel* model);
		virtual ~MyViscosity(void);

		static NonPressureForceBase* creator(FluidModel* model) { return new MyViscosity(model); }

		virtual void step();
		virtual void reset();

		virtual void performNeighborhoodSearchSort();

		FORCE_INLINE const Vector3r& getVDiff(const unsigned int i) const
		{
			return m_vDiff[i];
		}

		FORCE_INLINE Vector3r& getVDiff(const unsigned int i)
		{
			return m_vDiff[i];
		}

		FORCE_INLINE void setVDiff(const unsigned int i, const Vector3r& val)
		{
			m_vDiff[i] = val;
		}
	};
}

#endif