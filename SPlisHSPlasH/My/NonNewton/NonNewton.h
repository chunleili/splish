#pragma once

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
// #include "SPlisHSPlasH/NonPressureForceBase.h"
#include "SPlisHSPlasH/SurfaceTension/SurfaceTensionBase.h"

namespace SPH
{
	/** \brief NonNewton
	 * this class is based on XSPH
	 * 这个类的作用是实现非牛顿粘性。根据m_nonNewtonMethod的值来选择不同的非牛顿粘性模型。
	 * m_nonNewtonMethod = 0 时，使用牛顿粘性模型
	 * m_nonNewtonMethod = 1 时，使用Power law粘性模型
	 * m_nonNewtonMethod = 2 时，使用Cross模型
	 * m_nonNewtonMethod = 3 时，使用Casson模型
	*/
	class NonNewton : public SurfaceTensionBase
	{
	private:

		virtual void initParameters() override;
		void computeNonNewtonViscosity();
		void computeViscosityNewtonian();
		void computeViscosityPowerLaw();
		void computeViscosityCrossModel();
		void computeViscosityCassonModel();
		void computeViscosityCarreau();
		void computeViscosityBingham();
		void computeViscosityHerschelBulkley();
		void calcStrainRate();
		Real FNorm(const Vector6r& vec) const;
	public:
	    enum NonNewtonMethod { ENUM_NEWTONIAN=0, ENUM_POWER_LAW, ENUM_CROSS_MODEL, ENUM_CASSON_MODEL, ENUM_CARREAU, ENUM_BINGHAM, ENUM_HERSCHEL_BULKLEY};
		NonNewtonMethod m_nonNewtonMethod = ENUM_NEWTONIAN;

		std::vector<Vector6r> m_strainRate;
		std::vector<float> m_strainRateFNorm;
		std::vector<float> m_boundaryViscosity;
		std::vector<float> m_nonNewtonViscosity;

		float power_index = 0.3f;
		float consistency_index = 0.5;
		float m_viscosity0 = 10.0f; //initial viscosity
		float m_viscosity_inf = 1.0f; //infinite viscosity(must lower than initial viscosity)
		float m_muC = 0.00298f;
		float m_tauC = 0.02876f;
		float m_lambda = 4.020f;
		float m_maxViscosity = 0.0f;
		float m_avgViscosity = 0.0f;
		float critical_strainRate = 1e-2f;

		static int NON_NEWTON_METHOD;
		static int VISCOSITY_COEFFICIENT;
		static int VISCOSITY_COEFFICIENT_BOUNDARY;
		static int NEWTONIAN_;
		static int POWER_LAW_;
		static int CROSS_MODEL_;
		static int CASSON_MODEL_;
		static int POWER_INDEX;
		static int CONSISTENCY_INDEX;
		static int VISCOSITY0;
		static int VISCOSITY_INF;
		static int MU_C;
		static int TAU_C;
		static int LAMBDA;
		static int MAX_VISCOSITY;
		static int AVG_VISCOSITY;

		NonNewton(FluidModel *model);
        virtual void init() override;
		virtual ~NonNewton(void);
        virtual void step() override final;
        virtual void reset() override final;
        static NonPressureForceBase* creator(FluidModel* model) {return new NonNewton(model);}

		std::vector<float>& getViscosity() { return m_nonNewtonViscosity; }

		unsigned int numParticles;
	};
}