#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "SPlisHSPlasH/BoundaryModel_Koschier2017.h"
#include "SPlisHSPlasH/BoundaryModel_Bender2019.h"
#include "NonNewton.h"
#include "Utils/myPrint.h"
#include "Utils/mathUtils.h"

using namespace SPH;
using namespace GenParam;

NonNewton::NonNewton(FluidModel *model) :
SurfaceTensionBase(model)
{
	std::cout<<"constructor\n";

	numParticles = model->numActiveParticles();

	m_strainRateNorm.resize(numParticles, 0.0);
	m_nonNewtonViscosity.resize(numParticles, 0.0);
	m_strainRate.resize(numParticles, Vector6r::Zero());
	m_boundaryViscosity.resize(numParticles, 0.0);

	model->addField({ "strainRate", FieldType::Vector6, [&](const unsigned int i) -> Vector6r* { return &m_strainRate[i]; }, true });
	//nonNewtonViscosity already exists in FluidModel
	model->addField({ "strainRateNorm", FieldType::Scalar, [&](const unsigned int i) -> Real* { return &m_strainRateNorm[i]; }, true });

}

NonNewton::~NonNewton(void)
{
	m_model->removeFieldByName("strainRate");
	m_strainRate.clear();

	m_model->removeFieldByName("strainRateNorm");
	m_strainRateNorm.clear();
}

void NonNewton::init()
{
	std::cout<<"init\n";
	initParameters();
}

void NonNewton::initParameters()
{
	MAX_VISCOSITY = createNumericParameter("max_viscosity", "max_viscosity", &m_maxViscosity);
	AVG_VISCOSITY = createNumericParameter("average_viscosity", "average_viscosity", &m_avgViscosity);
	MIN_VISCOSITY = createNumericParameter("min_viscosity", "min_viscosity", &m_minViscosity);

	setGroup(MAX_VISCOSITY, "Viscosity");
	setGroup(AVG_VISCOSITY, "Viscosity");
	setGroup(MIN_VISCOSITY, "Viscosity");

	setDescription(MAX_VISCOSITY, "Max viscosity of all fluid particles.");
	setDescription(AVG_VISCOSITY, "Average viscosity of all fluid particles.");
	setDescription(MIN_VISCOSITY, "Min viscosity of all fluid particles.");

	getParameter(MAX_VISCOSITY)->setReadOnly(true);
	getParameter(AVG_VISCOSITY)->setReadOnly(true);
	getParameter(MIN_VISCOSITY)->setReadOnly(true);

	NON_NEWTON_METHOD = createEnumParameter("nonNewtonMethod", "nonNewtonMethod", &(int)m_nonNewtonMethod);
	setGroup(NON_NEWTON_METHOD, "Viscosity");
	setDescription(NON_NEWTON_METHOD, "Method for nonNewton.");
	EnumParameter *enumParam = static_cast<EnumParameter*>(getParameter(NON_NEWTON_METHOD));
	enumParam->addEnumValue("Newtonian", NEWTONIAN_);
	enumParam->addEnumValue("Power Law", POWER_LAW_);
	enumParam->addEnumValue("Cross Model", CROSS_);
	enumParam->addEnumValue("Casson Model", CASSON_);
	enumParam->addEnumValue("Carreau Model", CARREAU_);
	enumParam->addEnumValue("Bingham Model", BINGHAM_);
	enumParam->addEnumValue("Herschel-Bulkley Model", HERSCHEL_BULKLEY_);


	POWER_INDEX = createNumericParameter("power_index", "power_index", &power_index);
	CONSISTENCY_INDEX = createNumericParameter("consistency_index", "consistency_index", &consistency_index);
	VISCOSITY0 = createNumericParameter("viscosity0", "viscosity0", &m_viscosity0);
	VISCOSITY_INF = createNumericParameter("viscosity_inf", "viscosity_inf", &m_viscosity_inf);
	MU_C = createNumericParameter("muC", "muC", &m_muC);
	CRITICAL_STRAIN_RATE = createNumericParameter("criticalStrainRate", "criticalStrainRate", &m_criticalStrainRate);

	setGroup(POWER_INDEX, "Viscosity");
	setGroup(CONSISTENCY_INDEX, "Viscosity");
	setGroup(VISCOSITY0, "Viscosity");
	setGroup(VISCOSITY_INF, "Viscosity");
	setGroup(MU_C, "Viscosity");
	setGroup(CRITICAL_STRAIN_RATE, "Viscosity");

	setDescription(POWER_INDEX, "Power index for power law.");
	setDescription(CONSISTENCY_INDEX, "Consistency index for power law.");
	setDescription(VISCOSITY0, "Initial viscosity.");
	setDescription(VISCOSITY_INF, "Infinite viscosity for the cross model.");
	setDescription(MU_C, "Critical shear rate for the Casson model.");
	setDescription(CRITICAL_STRAIN_RATE, "Critical strain rate for the Herschel-Bulkley model.");

}



void NonNewton::reset()
{
	std::cout<<"reset!\n";

	numParticles = m_model->numActiveParticles();
	
	m_strainRate.resize(numParticles, Vector6r::Zero());
	m_nonNewtonViscosity.resize(numParticles, 0.0);
	m_boundaryViscosity.resize(numParticles, 0.0);
}

// void NonNewton::calcStrainRate()
// {
// // shear strain rate calculation
// 	Simulation *sim = Simulation::getCurrent();
// 	const unsigned int fluidModelIndex = m_model->getPointSetIndex();

// 	for (unsigned int i = 0; i < numParticles; ++i)
// 	{
// 		Vector6r strainRate;
// 		const Vector3r &xi = m_model->getPosition(i);
// 		const Vector3r &vi = m_model->getVelocity(i);
// 		const Real density_i = m_model->getDensity(i);

// 		for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i); j++)
// 		{
// 			const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, fluidModelIndex, i, j);
// 			const Vector3r &xj = m_model->getPosition(neighborIndex);

// 			const Vector3r &vj = m_model->getVelocity(neighborIndex);

// 			const Vector3r gradW = sim->gradW(xi - xj);
// 			const Vector3r vji = vj - vi;
// 			const Real m = m_model->getMass(neighborIndex);
// 			const Real m2 = m * static_cast<Real>(2.0);
// 			strainRate[0] += m2 * vji[0] * gradW[0];
// 			strainRate[1] += m2 * vji[1] * gradW[1];
// 			strainRate[2] += m2 * vji[2] * gradW[2];
// 			strainRate[3] += m * (vji[0] * gradW[1] + vji[1] * gradW[0]);
// 			strainRate[4] += m * (vji[0] * gradW[2] + vji[2] * gradW[0]);
// 			strainRate[5] += m * (vji[1] * gradW[2] + vji[2] * gradW[1]);
// 		}
// 		strainRate = (static_cast<Real>(0.5) / density_i) * strainRate;

// 		m_strainRate[i] = strainRate;

// 		m_strainRateNorm[i] = FNorm(m_strainRate[i]);
// 		// end shear strain rate calculation
// 	}
// }


void NonNewton::calcStrainRate()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();

	for (unsigned int i = 0; i < numParticles; ++i)
	{
		const Vector3r &xi = m_model->getPosition(i);
		const Vector3r &vi = m_model->getVelocity(i);
		const Real density_i = m_model->getDensity(i);

		Matrix3r velGrad = Matrix3r::Zero();
		Matrix3r strainRate = Matrix3r::Zero();

		for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i); j++)
		{
			const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, fluidModelIndex, i, j);
			const Vector3r &xj = m_model->getPosition(neighborIndex);
			const Vector3r &vj = m_model->getVelocity(neighborIndex);
			const Vector3r gradW = sim->gradW(xi - xj);
			const Vector3r vji = vj - vi;

			const Real m = m_model->getMass(neighborIndex);
			velGrad = vji * gradW.transpose();
			strainRate = velGrad + velGrad.transpose();
			strainRate *= m;
		}
		strainRate = (static_cast<Real>(0.5) / density_i) * strainRate;

		Real norm = 0.0f;
		norm = strainRate.norm();
		m_strainRateNorm[i] = norm;
	}
}

Real NonNewton::FNorm(const Vector6r & vec) const
{
	Real res = 0.0f;
	for (int i = 0; i < 6; i++)
		res += vec[i] * vec[i];
	res += vec[3] * vec[3];
	res += vec[4] * vec[4];
	res += vec[5] * vec[5];
	res = sqrt(res)	;
	return res;
}





void NonNewton::step()
{
	static int steps{0};
	steps++;
	// std::cout<<"\nstep: "<<steps<<"\n";
	numParticles = m_model->numActiveParticles();

	computeNonNewtonViscosity();

	//通过set函数将计算得到的粘度传递给fluidmodel
	for (unsigned int i = 0; i < numParticles; ++i)
	{
		m_model->setNonNewtonViscosity(i,m_nonNewtonViscosity[i]);
	}

	// [m_minViscosity, m_maxViscosity] = minMaxField(m_nonNewtonViscosity, numParticles);
	m_maxViscosity = maxField(m_nonNewtonViscosity, numParticles);
	m_minViscosity = minField(m_nonNewtonViscosity, numParticles);
	m_avgViscosity = avgField(m_nonNewtonViscosity, numParticles);
	// [m_minViscosity, m_maxViscosity, m_avgViscosity] = minMaxAvgField(m_boundaryViscosity, numParticles);

}

void NonNewton::computeNonNewtonViscosity()
{
	calcStrainRate();

	if(m_nonNewtonMethod == NonNewtonMethod::ENUM_POWER_LAW)
	{
		computeViscosityPowerLaw();
	}
	else if(m_nonNewtonMethod == NonNewtonMethod::ENUM_CROSS_MODEL)
	{
		computeViscosityCrossModel();
	}
	else if (m_nonNewtonMethod == NonNewtonMethod::ENUM_CASSON_MODEL)
	{
		computeViscosityCassonModel();
	}
	else if (m_nonNewtonMethod == NonNewtonMethod::ENUM_CARREAU)
	{
		computeViscosityCarreau();
	}
	else if (m_nonNewtonMethod == NonNewtonMethod::ENUM_BINGHAM)
	{
		computeViscosityBingham();
	}
	else if (m_nonNewtonMethod == NonNewtonMethod::ENUM_HERSCHEL_BULKLEY)
	{
		computeViscosityHerschelBulkley();
	}
	else
	{
		computeViscosityNewtonian();
	}

}



void NonNewton::computeViscosityNewtonian() 
{
	std::cout<<"computeViscosityNewtonian!\n";
	for (unsigned int i = 0; i < numParticles; ++i)
	{
		m_nonNewtonViscosity[i] = m_viscosity0;
	}
}

void NonNewton::computeViscosityPowerLaw() 
{
	for (unsigned int i = 0; i < numParticles; ++i)
	{
		m_nonNewtonViscosity[i] = consistency_index * pow(m_strainRateNorm[i], power_index - 1);
	}
}

void NonNewton::computeViscosityCrossModel() 
{
	assert((m_viscosity0 - m_viscosity_inf >= 0.0) && "the viscosity0 must be larger than viscosity_inf");
	if(m_viscosity0 - m_viscosity_inf < 0.0)
	{
		std::cout << "the viscosity0 must be larger than viscosity_inf" << std::endl;
		throw std::runtime_error("the viscosity0 must be larger than viscosity_inf");
	}
	for (unsigned int i = 0; i < numParticles; ++i)
	{

		m_nonNewtonViscosity[i] = m_viscosity_inf +  (m_viscosity0 - m_viscosity_inf) / (1 +  pow(consistency_index * m_strainRateNorm[i], power_index)) ;
	}
}

void NonNewton::computeViscosityCassonModel() 
{
	for (unsigned int i = 0; i < numParticles; ++i)
	{
		float res = sqrt(m_muC) +  sqrt(m_yieldStress / m_strainRateNorm[i]);
		m_nonNewtonViscosity[i] = res * res;
	}
}


void NonNewton::computeViscosityCarreau() 
{
	for (unsigned int i = 0; i < numParticles; ++i)
	{
		m_nonNewtonViscosity[i] = m_viscosity_inf +  (m_viscosity0 - m_viscosity_inf) / (1.0f +  pow(consistency_index * m_strainRateNorm[i]*m_strainRateNorm[i], (1.0f - power_index)/2.0f)) ;
	}
}



void NonNewton::computeViscosityBingham() 
{
	for (unsigned int i = 0; i < numParticles; ++i)
	{
		if(m_strainRateNorm[i] < m_criticalStrainRate)
			m_nonNewtonViscosity[i] = m_viscosity0;
		else
		{
			float tau0 = m_criticalStrainRate * (m_viscosity0 - m_viscosity_inf);
			m_nonNewtonViscosity[i] = tau0 / m_strainRateNorm[i] + m_viscosity_inf;
		}
	}
}

void NonNewton::computeViscosityHerschelBulkley() 
{
	for (unsigned int i = 0; i < numParticles; ++i)
	{
		if(m_strainRateNorm[i] < m_criticalStrainRate)
			m_nonNewtonViscosity[i] = m_viscosity0;
		else
		{
			float tau0 = m_viscosity0 * m_criticalStrainRate - consistency_index * pow(m_criticalStrainRate, power_index);
			m_nonNewtonViscosity[i] = tau0 / m_strainRateNorm[i] + consistency_index * pow(m_strainRateNorm[i], power_index - 1);
		}
	}
}


// void NonNewton::performNeighborhoodSearchSort()
// {
// 	Simulation *sim = Simulation::getCurrent();
// 	auto const& d = sim->getNeighborhoodSearch()->point_set(m_model->getPointSetIndex());
// 	d.sort_field(&m_nonNewtonViscosity[0]);
// } 