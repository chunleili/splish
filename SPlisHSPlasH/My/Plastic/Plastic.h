#ifndef __Plastic_h__
#define __Plastic_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SPlisHSPlasH/Elasticity/ElasticityBase.h"

namespace SPH
{
	/** \brief 
	 * Implemetation of plasticity
	 * Modified from the Elatticity_Becker2009
	*/
	class Plastic : public ElasticityBase
	{
	protected:
		// initial particle indices, used to access their original positions
		std::vector<unsigned int> m_current_to_initial_index;
		std::vector<unsigned int> m_initial_to_current_index;
		// initial particle neighborhood
		std::vector<std::vector<unsigned int>> m_initialNeighbors;
		// volumes in rest configuration
		std::vector<Real> m_restVolumes;
		std::vector<Matrix3r> m_rotations;
		std::vector<Vector6r> m_stress;
		std::vector<Matrix3r> m_F;
		Real m_alpha;

		std::vector<Vector6r> m_plasticStrain; //add Plastic Strain 

		void initValues();
		void computeRotations();
		void computeStress();
		void computeForces();

		void computePlasticStrain(Vector6r & totalStrain);
		void computeNablaU(int i, Matrix3r &nablaU);
		void computeTotalStrain(Matrix3r &nablaU, Vector6r & totalStrain);

		virtual void initParameters();

		//////////////////////////////////////////////////////////////////////////
		// multiplication of symmetric matrix, represented by a 6D vector, and a 
		// 3D vector
		//////////////////////////////////////////////////////////////////////////
		FORCE_INLINE void symMatTimesVec(const Vector6r & M, const Vector3r & v, Vector3r &res)
		{
			res[0] = M[0] * v[0] + M[3] * v[1] + M[4] * v[2];
			res[1] = M[3] * v[0] + M[1] * v[1] + M[5] * v[2];
			res[2] = M[4] * v[0] + M[5] * v[1] + M[2] * v[2];
		}


	public:
		static int ALPHA;

		//user defined parameters
		float elasticLimit;
		float plasticLimit;
		static int ELASTIC_LIMIT;
		static int PLASTIC_LIMIT;

		Plastic(FluidModel *model);
		virtual ~Plastic(void);

		static NonPressureForceBase* creator(FluidModel* model) { return new Plastic(model); }

		virtual void step();
		virtual void reset();
		virtual void performNeighborhoodSearchSort();

		virtual void saveState(BinaryFileWriter &binWriter);
		virtual void loadState(BinaryFileReader &binReader);
	};
}

#endif
