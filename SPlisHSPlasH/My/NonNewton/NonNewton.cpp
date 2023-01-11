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

int NonNewton::VISCOSITY_COEFFICIENT_BOUNDARY = -1;
int NonNewton::VISCOSITY_COEFFICIENT = -1;
int NonNewton::NON_NEWTON_METHOD = -1;
int NonNewton::NEWTONIAN_ = -1;
int NonNewton::POWER_LAW_ = -1;
int NonNewton::CROSS_MODEL_ = -1;
int NonNewton::CASSON_MODEL_ = -1;
int NonNewton::POWER_INDEX = -1;
int NonNewton::CONSISTENCY_INDEX = -1;
int NonNewton::VISCOSITY0 = -1;
int NonNewton::VISCOSITY_INF = -1;
int NonNewton::MU_C = -1;
int NonNewton::TAU_C = -1;
int NonNewton::LAMBDA = -1;
int NonNewton::MAX_VISCOSITY = -1;
int NonNewton::AVG_VISCOSITY = -1;


NonNewton::NonNewton(FluidModel *model) :
SurfaceTensionBase(model)
{
	std::cout<<"constructor\n";
	// m_nonNewtonMethod = 0;
	// power_index = 0.5;
	// consistency_index = 1.0;
	// m_viscosity0 = 10.0f;
	// m_viscosity_inf = 1.0f;
	// m_muC = 0.00298;
	// m_tauC = 0.02876;
	// m_lambda = 4.020;

	numParticles = model->numActiveParticles();

	m_strainRateFNorm.resize(numParticles, 0.0);

	m_strainRate.resize(numParticles, Vector6r::Zero());
	model->addField({ "strainRate", FieldType::Vector6, [&](const unsigned int i) -> Vector6r* { return &m_strainRate[i]; }, true });

	m_nonNewtonViscosity.resize(numParticles, 0.0);
	model->addField({ "viscosity", FieldType::Scalar, [&](const unsigned int i) -> Real* { return &m_nonNewtonViscosity[i]; }, true });

	m_boundaryViscosity.resize(numParticles, 0.0);
}

NonNewton::~NonNewton(void)
{
	m_model->removeFieldByName("strainRate");
	m_strainRate.clear();

	m_model->removeFieldByName("viscosity");
	m_nonNewtonViscosity.clear();
}

void NonNewton::init()
{
	std::cout<<"init\n";
	initParameters();
}

void NonNewton::initParameters()
{
	MAX_VISCOSITY = createNumericParameter("max_viscosity", "max_viscosity", &m_maxViscosity);
	setGroup(MAX_VISCOSITY, "Viscosity");
	setDescription(MAX_VISCOSITY, "Max viscosity of all fluid particles.");
	getParameter(MAX_VISCOSITY)->setReadOnly(true);

	AVG_VISCOSITY = createNumericParameter("average_viscosity", "average_viscosity", &m_avgViscosity);
	setGroup(AVG_VISCOSITY, "Viscosity");
	setDescription(AVG_VISCOSITY, "Average viscosity of all fluid particles.");
	getParameter(AVG_VISCOSITY)->setReadOnly(true);

	NON_NEWTON_METHOD = createEnumParameter("nonNewtonMethod", "nonNewtonMethod", &(int)m_nonNewtonMethod);
	setGroup(NON_NEWTON_METHOD, "Viscosity");
	setDescription(NON_NEWTON_METHOD, "Method for nonNewton.");
	EnumParameter *enumParam = static_cast<EnumParameter*>(getParameter(NON_NEWTON_METHOD));
	enumParam->addEnumValue("Newtonian", NEWTONIAN_);
	enumParam->addEnumValue("Power Law", POWER_LAW_);
	enumParam->addEnumValue("Cross Model", CROSS_MODEL_);
	enumParam->addEnumValue("Casson Model", CASSON_MODEL_);

	POWER_INDEX = createNumericParameter("power_index", "power_index", &power_index);
	CONSISTENCY_INDEX = createNumericParameter("consistency_index", "consistency_index", &consistency_index);
	setGroup(POWER_INDEX, "Viscosity");
	setGroup(CONSISTENCY_INDEX, "Viscosity");
	setDescription(POWER_INDEX, "Power index for power law.");
	setDescription(CONSISTENCY_INDEX, "Consistency index for power law.");

	VISCOSITY0 = createNumericParameter("viscosity0", "viscosity0", &m_viscosity0);
	VISCOSITY_INF = createNumericParameter("viscosity_inf", "viscosity_inf", &m_viscosity_inf);
	setGroup(VISCOSITY0, "Viscosity");
	setGroup(VISCOSITY_INF, "Viscosity");
	setDescription(VISCOSITY0, "Initial viscosity.");
	setDescription(VISCOSITY_INF, "Infinite viscosity for the cross model.");

	MU_C = createNumericParameter("muC", "muC", &m_muC);
	TAU_C = createNumericParameter("tauC", "tauC", &m_tauC);
	LAMBDA = createNumericParameter("lambda", "lambda", &m_lambda);
	setGroup(MU_C, "Viscosity");
	setGroup(TAU_C, "Viscosity");
	setGroup(LAMBDA, "Viscosity");
}

void NonNewton::computeNonNewtonViscosity()
{
	if(m_nonNewtonMethod == NonNewtonMethod::ENUM_POWER_LAW)
	{
		calcStrainRate();
		computeViscosityPowerLaw();
	}
	else if(m_nonNewtonMethod == NonNewtonMethod::ENUM_CROSS_MODEL)
	{
		calcStrainRate();
		computeViscosityCrossModel();
	}
	else if (m_nonNewtonMethod == NonNewtonMethod::ENUM_CASSON_MODEL)
	{
		calcStrainRate();
		computeViscosityCassonModel();
	}
	else if (m_nonNewtonMethod == NonNewtonMethod::ENUM_CARREAU)
	{
		calcStrainRate();
		computeViscosityCarreau();
	}
	else if (m_nonNewtonMethod == NonNewtonMethod::ENUM_BINGHAM)
	{
		calcStrainRate();
		computeViscosityBingham();
	}
	else if (m_nonNewtonMethod == NonNewtonMethod::ENUM_HERSCHEL_BULKLEY)
	{
		calcStrainRate();
		computeViscosityHerschelBulkley();
	}
	else
	{
		// std::cout<<"Newtonian!\n";
		computeViscosityNewtonian();
	}

	m_maxViscosity = 0.0;
	m_maxViscosity = maxField(m_nonNewtonViscosity, m_model->numActiveParticles());
	echo(m_maxViscosity);

	m_avgViscosity = averageField(m_nonNewtonViscosity, m_model->numActiveParticles());
	echo(m_avgViscosity);
}



void NonNewton::step()
{
	static int steps{0};
	steps++;
	std::cout<<"\nstep: "<<steps<<"\n";
	numParticles = m_model->numActiveParticles();

	computeNonNewtonViscosity();

	//通过set函数将计算得到的粘度传递给fluidmodel
	for (unsigned int i = 0; i < numParticles; ++i)
	{
		m_model->setNonNewtonViscosity(i,m_nonNewtonViscosity[i]);
	}
}




void NonNewton::reset()
{
	std::cout<<"reset!\n";

	// m_nonNewtonMethod = 0;
	// power_index = 0.5;
	// consistency_index = 1.0;
	// m_viscosity0 = 10.0f;
	// m_viscosity_inf = 1.0f;
	// m_muC = 0.00298;
	// m_tauC = 0.02876;
	// m_lambda = 4.020;
	// m_maxViscosity = 0.0;

	numParticles = m_model->numActiveParticles();
	
	m_strainRate.resize(numParticles, Vector6r::Zero());
	m_nonNewtonViscosity.resize(numParticles, 0.0);
	m_boundaryViscosity.resize(numParticles, 0.0);
}

void NonNewton::calcStrainRate()
{
// shear strain rate calculation
	Simulation *sim = Simulation::getCurrent();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();

	for (unsigned int i = 0; i < numParticles; ++i)
	{
		Vector6r strainRate;
		const Vector3r &xi = m_model->getPosition(i);
		const Vector3r &vi = m_model->getVelocity(i);
		const Real density_i = m_model->getDensity(i);

		for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i); j++)
		{
			const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, fluidModelIndex, i, j);
			const Vector3r &xj = m_model->getPosition(neighborIndex);

			const Vector3r &vj = m_model->getVelocity(neighborIndex);

			const Vector3r gradW = sim->gradW(xi - xj);
			const Vector3r vji = vj - vi;
			const Real m = m_model->getMass(neighborIndex);
			const Real m2 = m * static_cast<Real>(2.0);
			strainRate[0] += m2 * vji[0] * gradW[0];
			strainRate[1] += m2 * vji[1] * gradW[1];
			strainRate[2] += m2 * vji[2] * gradW[2];
			strainRate[3] += m * (vji[0] * gradW[1] + vji[1] * gradW[0]);
			strainRate[4] += m * (vji[0] * gradW[2] + vji[2] * gradW[0]);
			strainRate[5] += m * (vji[1] * gradW[2] + vji[2] * gradW[1]);
		}
		strainRate = (static_cast<Real>(0.5) / density_i) * strainRate;

		m_strainRate[i] = strainRate;

		m_strainRateFNorm[i] = FNorm(m_strainRate[i]);
		// end shear strain rate calculation
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

void NonNewton::computeViscosityNewtonian() 
{
	for (unsigned int i = 0; i < numParticles; ++i)
	{
		m_nonNewtonViscosity[i] = m_viscosity0;
	}
}

void NonNewton::computeViscosityPowerLaw() 
{
	
	for (unsigned int i = 0; i < numParticles; ++i)
	{
		m_nonNewtonViscosity[i] = consistency_index * pow(m_strainRateFNorm[i], power_index - 1);
	}
}

void NonNewton::computeViscosityCrossModel() 
{
	
	assert(m_viscosity0 - m_viscosity_inf >= 0.0);
	if(m_viscosity0 - m_viscosity_inf < 0.0)
	{
		std::cout << "the viscosity0 must be larger than viscosity_inf" << std::endl;
		throw std::runtime_error("the viscosity0 must be larger than viscosity_inf");
	}

	// echo(m_strainRate[0].norm());
	// std::cout<<"pow(m_strainRate[0].norm(), power_index)"<<pow(m_strainRate[0].norm(), power_index)<<std::endl;

	for (unsigned int i = 0; i < numParticles; ++i)
	{

		m_nonNewtonViscosity[i] = m_viscosity_inf +  (m_viscosity0 - m_viscosity_inf) / (1 +  pow(consistency_index * m_strainRateFNorm[i], power_index)) ;
	}
}

void NonNewton::computeViscosityCassonModel() 
{
	
	for (unsigned int i = 0; i < numParticles; ++i)
	{
		m_nonNewtonViscosity[i] = sqrt(m_muC) + sqrt(m_tauC) / (sqrt(m_lambda) + sqrt(m_strainRateFNorm[i]));
	}
}


void NonNewton::computeViscosityCarreau() 
{
	
	for (unsigned int i = 0; i < numParticles; ++i)
	{
		m_nonNewtonViscosity[i] = m_viscosity_inf +  (m_viscosity0 - m_viscosity_inf) / (1.0f +  pow(consistency_index * FNorm(m_strainRate[i])*FNorm(m_strainRate[i]), (1.0f - power_index)/2.0)) ;
	}
}

void NonNewton::computeViscosityBingham() 
{
	
	for (unsigned int i = 0; i < numParticles; ++i)
	{
		if(FNorm(m_strainRate[i]) < critical_strainRate)
			m_nonNewtonViscosity[i] = m_viscosity0;
		else
		{
			float tau0 = critical_strainRate * (m_viscosity0 - m_viscosity_inf);
			m_nonNewtonViscosity[i] = tau0 / FNorm(m_strainRate[i]) + m_viscosity_inf;
		}
	}
}

void NonNewton::computeViscosityHerschelBulkley() 
{
	
	for (unsigned int i = 0; i < numParticles; ++i)
	{
		if(FNorm(m_strainRate[i]) < critical_strainRate)
			m_nonNewtonViscosity[i] = m_viscosity0;
		else
		{
			float tau0 = critical_strainRate * (m_viscosity0 - m_viscosity_inf);
			m_nonNewtonViscosity[i] = tau0 / FNorm(m_strainRate[i]) + consistency_index * pow(FNorm(m_strainRate[i]), power_index - 1);
		}
	}
}