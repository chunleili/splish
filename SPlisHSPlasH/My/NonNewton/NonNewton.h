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
	 * m_nonNewtonMethod = 4 时，使用Carreau模型
	 * m_nonNewtonMethod = 5 时，使用Bingham模型
	 * m_nonNewtonMethod = 6 时，使用Herschel-Bulkley模型
	 * 
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
		std::vector<float> m_strainRateNorm;
		std::vector<float> m_boundaryViscosity;
		std::vector<float> m_nonNewtonViscosity;

		float power_index = 0.667f;
		float consistency_index = 100.0f;
		float m_viscosity0 = 2000.0f; //initial viscosity
		float m_viscosity_inf = 1.0f; //infinite viscosity(must lower than initial viscosity)
		float m_criticalStrainRate = 20.0f;
		float m_muC = 10.0f;
		// float m_yieldStress = (m_viscosity0 - m_viscosity_inf) * m_criticalStrainRate;
		float m_yieldStress = 200.0f;
		// float m_lambda = 0.5;
		float m_maxViscosity = 0.0f;
		float m_avgViscosity = 0.0f;
		float m_minViscosity = 0.0f;

		inline static int VISCOSITY_COEFFICIENT;
		inline static int VISCOSITY_COEFFICIENT_BOUNDARY;

		inline static int POWER_INDEX = -1;
		inline static int CONSISTENCY_INDEX = -1;
		inline static int VISCOSITY0 = -1;
		inline static int VISCOSITY_INF = -1;
		inline static int MU_C = -1;

		inline static int MAX_VISCOSITY = -1;
		inline static int AVG_VISCOSITY = -1;
		inline static int MIN_VISCOSITY = -1;

		inline static int NON_NEWTON_METHOD;
		inline static int NEWTONIAN_ = -1;
		inline static int POWER_LAW_ = -1;
		inline static int CROSS_ = -1;
		inline static int CASSON_	= -1;
		inline static int CARREAU_ = -1;
		inline static int BINGHAM_ = -1;
		inline static int HERSCHEL_BULKLEY_ = -1;


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