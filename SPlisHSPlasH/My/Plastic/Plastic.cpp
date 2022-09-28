#include "Plastic.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/Utilities/MathFunctions.h"
#include "../Utils/myPrint.h"

using namespace SPH;
using namespace GenParam;

int Plastic::ALPHA = -1;
int Plastic::ELASTIC_LIMIT = -1;
int Plastic::PLASTIC_LIMIT = -1;


Plastic::Plastic(FluidModel *model) :
	ElasticityBase(model)
{
	numParticles = model->numActiveParticles();
	m_restVolumes.resize(numParticles);
	m_current_to_initial_index.resize(numParticles);
	m_initial_to_current_index.resize(numParticles);
	m_initialNeighbors.resize(numParticles);
	m_rotations.resize(numParticles, Matrix3r::Identity());
	m_stress.resize(numParticles);
	m_F.resize(numParticles);
	m_alpha = 0.0;

	elasticLimit = 0.001;
	plasticLimit = 0.486;

	m_plasticStrain.resize(numParticles);
	m_isFracture.resize(numParticles,0);
	m_totalStrain.resize(numParticles);
	m_elasticStrain.resize(numParticles);

	initValues();

	model->addField({ "rest volume", FieldType::Scalar, [&](const unsigned int i) -> Real* { return &m_restVolumes[i]; }, true });
	model->addField({ "rotation", FieldType::Matrix3, [&](const unsigned int i) -> Real* { return &m_rotations[i](0,0); } });
	model->addField({ "stress", FieldType::Vector6, [&](const unsigned int i) -> Real* { return &m_stress[i][0]; } });
	model->addField({ "deformation gradient", FieldType::Matrix3, [&](const unsigned int i) -> Real* { return &m_F[i](0,0); } });

	model->addField({ "plastic strain", FieldType::Vector6, [&](const unsigned int i) -> Real* { return &m_plasticStrain[i][0]; }, true });
	model->addField({ "isFracture", FieldType::UInt, [&](const unsigned int i) -> int* { return &m_isFracture[i]; }, true  });
}

Plastic::~Plastic(void)
{
	m_model->removeFieldByName("rest volume");
	m_model->removeFieldByName("rotation");
	m_model->removeFieldByName("stress");
	m_model->removeFieldByName("deformation gradient");

	m_model->removeFieldByName("plastic strain");
	m_model->removeFieldByName("isFracture");

}

void Plastic::initParameters()
{
	ElasticityBase::initParameters();

	ALPHA = createNumericParameter("alpha", "Zero-energy modes suppression", &m_alpha);
	setGroup(ALPHA, "Elasticity");
	setDescription(ALPHA, "Coefficent for zero-energy modes suppression method");
	RealParameter *rparam = static_cast<RealParameter*>(getParameter(ALPHA));
	rparam->setMinValue(0.0);

	//for plasticity
	ELASTIC_LIMIT = createNumericParameter("elasticLimit", "elastic limit", &elasticLimit);
	setGroup(ELASTIC_LIMIT, "Plastic");
	setDescription(ELASTIC_LIMIT, "elastic limit");
	rparam = static_cast<RealParameter*>(getParameter(ELASTIC_LIMIT));
	rparam->setMinValue(0.0);

	PLASTIC_LIMIT = createNumericParameter("plasticLimit", "plastic limit", &plasticLimit);
	setGroup(PLASTIC_LIMIT, "Plastic");
	setDescription(PLASTIC_LIMIT, "plastic limit");
	rparam = static_cast<RealParameter*>(getParameter(PLASTIC_LIMIT));
	rparam->setMinValue(0.0);
}

void Plastic::initValues()
{
	Simulation *sim = Simulation::getCurrent();
	sim->getNeighborhoodSearch()->find_neighbors();

	FluidModel *model = m_model;
	const unsigned int fluidModelIndex = model->getPointSetIndex();

	// Store the neighbors in the reference configurations and
	// compute the volume of each particle in rest state
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static) 
		for (int i = 0; i < (int)numParticles; i++)
		{
			m_current_to_initial_index[i] = i;
			m_initial_to_current_index[i] = i;

			// only neighbors in same phase will influence elasticity
			const unsigned int numNeighbors = sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i);
			m_initialNeighbors[i].resize(numNeighbors);
			for (unsigned int j = 0; j < numNeighbors; j++)
				m_initialNeighbors[i][j] = sim->getNeighbor(fluidModelIndex, fluidModelIndex, i, j);

			// compute volume
			Real density = model->getMass(i) * sim->W_zero();
			const Vector3r &xi = model->getPosition(i);
			forall_fluid_neighbors_in_same_phase(
				density += model->getMass(neighborIndex) * sim->W(xi - xj);
			)
			m_restVolumes[i] = model->getMass(i) / density;

			//setzero for plastic strain
			m_plasticStrain[i].setZero();
			m_elasticStrain[i].setZero();
		}
	}
}

void Plastic::step()
{
	// computeRotations();
	computeStress();
	computeForces();
	m_step++;
	// std::string fname = "m_initial_to_current_index_" + std::to_string(m_step)+ ".txt";
	// printScalarField(fname,m_initial_to_current_index,0);
	// std::string fname = "m_initialNeighbors_" + std::to_string(m_step)+ ".txt";
	// printVectorField(fname,m_initialNeighbors,0);

}


void Plastic::reset()
{
	initValues();
}

void Plastic::performNeighborhoodSearchSort()
{
	if (numParticles == 0)
		return;

	Simulation *sim = Simulation::getCurrent();
	auto const& d = sim->getNeighborhoodSearch()->point_set(m_model->getPointSetIndex());
	d.sort_field(&m_restVolumes[0]);
	d.sort_field(&m_current_to_initial_index[0]);
	d.sort_field(&m_plasticStrain[0]);
	d.sort_field(&m_isFracture[0]);

	for (unsigned int i = 0; i < numParticles; i++)
		m_initial_to_current_index[m_current_to_initial_index[i]] = i;
}

void Plastic::computeRotations()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	FluidModel *model = m_model;

// #pragma omp parallel default(shared)
	{
// #pragma omp for schedule(static)
		for (int i = 0; i < (int)numParticles; i++)
		{
			if (model->getParticleState(i) == ParticleState::Active && 		
				m_isFracture[i] == 0)
			{

				const unsigned int i0 = m_current_to_initial_index[i];
				const Vector3r &xi = m_model->getPosition(i);
				const Vector3r &xi0 = m_model->getPosition0(i0);
				Matrix3r Apq;
				Apq.setZero();

				const size_t numNeighbors = m_initialNeighbors[i0].size();

				//////////////////////////////////////////////////////////////////////////
				// Fluid
				//////////////////////////////////////////////////////////////////////////
				for (unsigned int j = 0; j < numNeighbors; j++)
				{
					const unsigned int neighborIndex = m_initial_to_current_index[m_initialNeighbors[i0][j]];
					// get initial neighbor index considering the current particle order
					const unsigned int neighborIndex0 = m_initialNeighbors[i0][j];

					const Vector3r &xj = model->getPosition(neighborIndex);
					const Vector3r &xj0 = m_model->getPosition0(neighborIndex0);
					const Vector3r xj_xi = xj - xi;
					const Vector3r xj_xi_0 = xj0 - xi0;
					Apq += m_model->getMass(neighborIndex) * sim->W(xj_xi_0) * (xj_xi * xj_xi_0.transpose());
				}

				// 			Vector3r sigma;
				// 			Matrix3r U, VT;
				// 			MathFunctions::svdWithInversionHandling(Apq, sigma, U, VT);
				// 			m_rotations[i] = U * VT;
				Quaternionr q(m_rotations[i]);
				MathFunctions::extractRotation(Apq, q, 10);
				m_rotations[i] = q.matrix();
			}
		}
	}
}

void Plastic::computeStress()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	FluidModel *model = m_model;

	// Elasticity tensor
	Matrix6r C;
	C.setZero();
	const Real factor = m_youngsModulus / ((static_cast<Real>(1.0) + m_poissonRatio)*(static_cast<Real>(1.0) - static_cast<Real>(2.0) * m_poissonRatio));
	C(0, 0) = C(1, 1) = C(2, 2) = factor * (static_cast<Real>(1.0) - m_poissonRatio);
	C(0, 1) = C(0, 2) = C(1, 0) = C(1, 2) = C(2, 0) = C(2, 1) = factor * (m_poissonRatio);
	C(3, 3) = C(4, 4) = C(5, 5) = factor * static_cast<Real>(0.5)*(static_cast<Real>(1.0) - static_cast<Real>(2.0) * m_poissonRatio);

	// #pragma omp parallel default(shared)
	{
		// #pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			if (model->getParticleState(i) == ParticleState::Active
				&& m_isFracture[i] == 0
			)
			{
				Matrix3r nablaU;
				computeNablaU(i, nablaU);

				m_F[i] = nablaU + Matrix3r::Identity();
				
				Vector6r totalStrain;
				computeTotalStrain(nablaU, totalStrain);
				m_totalStrain[i] = totalStrain;

				computePlasticStrain(i, m_elasticStrain[i]); //FIXME: with bug

				Vector6r elasticStrain;
				elasticStrain = totalStrain - m_plasticStrain[i];
				m_elasticStrain[i] = elasticStrain;
				m_stress[i] = C * elasticStrain;
			}
			else
				m_stress[i].setZero();
		}
	}
}


/**
 * @brief compute NablaU
 * 
 * @param i: input, particle index
 * @param nablaU: output, gradient of displacement
 */
void Plastic::computeNablaU(int i, Matrix3r &nablaU)
{
	FluidModel *model = m_model;
	Simulation *sim = Simulation::getCurrent();

	const unsigned int i0 = m_current_to_initial_index[i];
	const Vector3r& xi = m_model->getPosition(i);
	const Vector3r& xi0 = m_model->getPosition0(i0);

	nablaU.setZero();
	const size_t numNeighbors = m_initialNeighbors[i0].size();

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	for (unsigned int j = 0; j < numNeighbors; j++)
	{
		const unsigned int neighborIndex = m_initial_to_current_index[m_initialNeighbors[i0][j]];
		// get initial neighbor index considering the current particle order 
		const unsigned int neighborIndex0 = m_initialNeighbors[i0][j];

		const Vector3r& xj = model->getPosition(neighborIndex);
		const Vector3r& xj0 = m_model->getPosition0(neighborIndex0);

		const Vector3r xj_xi = xj - xi;
		const Vector3r xj_xi_0 = xj0 - xi0;

		const Vector3r uji = m_rotations[i].transpose() * xj_xi - xj_xi_0;

		// m_displacement[i] = uji;
		// subtract because kernel gradient is taken in direction of xji0 instead of xij0
		// Eq(6) in Becker2009
		nablaU -= (m_restVolumes[neighborIndex] * uji) * sim->gradW(xj_xi_0).transpose();
	}	
}

/**
 * @brief compute totalStrain
 * 
 * @param nablaU: input, gradient of displacement
 * @param totalStrain: output
 */
void Plastic::computeTotalStrain(Matrix3r &nablaU, Vector6r & totalStrain)
{
	totalStrain.setZero();
	// compute Cauchy strain: epsilon = 0.5 (nabla u + nabla u^T)
	totalStrain[0] = nablaU(0, 0);						// \epsilon_{00}
	totalStrain[1] = nablaU(1, 1);						// \epsilon_{11}
	totalStrain[2] = nablaU(2, 2);						// \epsilon_{22}
	totalStrain[3] = static_cast<Real>(0.5) * (nablaU(0, 1) + nablaU(1, 0)); // \epsilon_{01}
	totalStrain[4] = static_cast<Real>(0.5) * (nablaU(0, 2) + nablaU(2, 0)); // \epsilon_{02}
	totalStrain[5] = static_cast<Real>(0.5) * (nablaU(1, 2) + nablaU(2, 1)); // \epsilon_{12}
}

/**
 * @brief Compute the plastic strain based on O'Brien 2002. 
 * 
 * Ref: James F. O’Brien et. al. 2002, "Graphical Modeling and Animation of Ductile Fracture"
 * 
 * @param i: input
 */
void Plastic::computePlasticStrain(int i, Vector6r &strain)
{
	//Eq(2) in O'Brien 2002
	Vector6r deviation =  strain;
	Real trace = strain[0] + strain[1] + strain[2];
	trace /= 3.0;
	deviation[0] -= trace;
	deviation[1] -= trace;
	deviation[2] -= trace;

	//Eq(3),  Consider plasticity only if exceeding elasticLimit
	Real deviationNorm = FNorm(deviation);
	// if(i==0)
		// printf("deviationNorm: %.4f, step %d\n", i, m_step);
	if(deviationNorm <= elasticLimit)
	{
		return;
	}
	// printf("particle %d is plastic in step %d\n", i, m_step);
	
	//Eq(4)
	Vector6r strainIncrement  = ( deviationNorm - elasticLimit ) / deviationNorm * deviation;

	//Eq(5) Calculate the accumulated plastic strain
	Vector6r newPlastic = m_plasticStrain[i] + strainIncrement;
	Real plasticNorm = FNorm(newPlastic);
	// if(deviationNorm >= plasticLimit) //barrier of plastic strain
	// {
	// 	m_isFracture[i] = 1;
	// 	printf("particle %d is fractured in step %d\n", i, m_step);
	// }
	Real ratio = (plasticLimit / plasticNorm);
	Real min = 1.0 < ratio ? 1.0 : ratio;
	m_plasticStrain[i] = newPlastic * min;
}



void Plastic::computeForces()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int numParticles = m_model->numActiveParticles();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	FluidModel *model = m_model;

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			if (model->getParticleState(i) == ParticleState::Active
				&& m_isFracture[i] == 0)
			{
				const unsigned int i0 = m_current_to_initial_index[i];
				const Vector3r& xi0 = m_model->getPosition0(i0);

				const size_t numNeighbors = m_initialNeighbors[i0].size();
				Vector3r fi;
				fi.setZero();

				//////////////////////////////////////////////////////////////////////////
				// Fluid
				//////////////////////////////////////////////////////////////////////////
				for (unsigned int j = 0; j < numNeighbors; j++)
				{
					const unsigned int neighborIndex = m_initial_to_current_index[m_initialNeighbors[i0][j]];
					// get initial neighbor index considering the current particle order 
					const unsigned int neighborIndex0 = m_initialNeighbors[i0][j];

					const Vector3r& xj0 = m_model->getPosition0(neighborIndex0);

					const Vector3r xj_xi_0 = xj0 - xi0;
					const Vector3r gradW0 = sim->gradW(xj_xi_0);

					const Vector3r dji = m_restVolumes[i] * gradW0;
					const Vector3r dij = -m_restVolumes[neighborIndex] * gradW0;

					Vector3r sdji, sdij;
					symMatTimesVec(m_stress[neighborIndex], dji, sdji); //stress times dji
					symMatTimesVec(m_stress[i], dij, sdij);

					const Vector3r fij = -m_restVolumes[neighborIndex] * sdji;
					const Vector3r fji = -m_restVolumes[i] * sdij;

					fi += m_rotations[neighborIndex] * fij - m_rotations[i] * fji;
				}
				fi = 0.5*fi;

				// elastic acceleration
				Vector3r& ai = model->getAcceleration(i);
				ai += fi / model->getMass(i);
			}
		}
	}
}

void SPH::Plastic::saveState(BinaryFileWriter &binWriter)
{
	binWriter.writeBuffer((char*)m_current_to_initial_index.data(), m_current_to_initial_index.size() * sizeof(unsigned int));
	binWriter.writeBuffer((char*)m_initial_to_current_index.data(), m_initial_to_current_index.size() * sizeof(unsigned int));
}

void SPH::Plastic::loadState(BinaryFileReader &binReader)
{
	binReader.readBuffer((char*)m_current_to_initial_index.data(), m_current_to_initial_index.size() * sizeof(unsigned int));
	binReader.readBuffer((char*)m_initial_to_current_index.data(), m_initial_to_current_index.size() * sizeof(unsigned int));
}

