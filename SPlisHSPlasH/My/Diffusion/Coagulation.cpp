#include "Coagulation.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/TimeManager.h"

using namespace SPH;
using namespace GenParam;

int Coagulation::THRES_HIGH  = -1;
int Coagulation::THRES_LOW = -1;
int Coagulation::DIFFUSIVITY = -1;
int Coagulation::R_SOURCE = -1;
int Coagulation::BOX_MIN = -1;
int Coagulation::BOX_MAX = -1;
int Coagulation::POINT_SRC_VAL = -1;
int Coagulation::POINT_SRC_POS = -1;

Coagulation::Coagulation(FluidModel* model) :
    NonPressureForceBase(model), 
    m_ccf(),
    m_thresHigh(static_cast<Real>(1.0)),
    m_thresLow(static_cast<Real>(0.1)),
    m_diffusivity(static_cast<Real>(50.0)),
    m_rSource(static_cast<Real>(0.1))
{
    m_ccf.resize(model->numParticles(), 0.0);

    model->addField({ "ccf field", FieldType::Scalar, [&](const unsigned int i) -> Real* { return &m_ccf[i]; } });

	m_boxMin.setZero();
	m_boxMax = Vector3r(0.3, 0.3, 0.3);

    // m_ccf[1501] = 100.0;
}



Coagulation::~Coagulation(void)
{
    m_model->removeFieldByName("ccf field");
    m_ccf.clear();
}

void Coagulation::initParameters()
{
    THRES_HIGH = createNumericParameter("thresHigh", "threshold high", &m_thresHigh);
	setGroup(THRES_HIGH, "coagualtion");
	setDescription(THRES_HIGH, "higher treshold for coagualtion");
	GenParam::RealParameter* rparam = static_cast<GenParam::RealParameter*>(getParameter(THRES_HIGH));
	rparam->setMinValue(0.0);

    THRES_LOW = createNumericParameter("thresLow", "threshold low", &m_thresLow);
    setGroup(THRES_LOW, "coagualtion");
    setDescription(THRES_LOW, "lower treshold for coagualtion");
    rparam = static_cast<GenParam::RealParameter*>(getParameter(THRES_LOW));
    rparam->setMinValue(0.0);

    DIFFUSIVITY = createNumericParameter("diffusivity", "diffusivity", &m_diffusivity);
    setGroup(DIFFUSIVITY, "coagualtion");
    setDescription(DIFFUSIVITY, "diffusivity of the convection-diffusion equation");
    rparam = static_cast<GenParam::RealParameter*>(getParameter(DIFFUSIVITY));
    rparam->setMinValue(0.0);

    R_SOURCE = createNumericParameter("rSource", "rSource", &m_rSource);
    setGroup(R_SOURCE, "coagualtion");
    setDescription(R_SOURCE, "source term of the convection-diffusion equation");
    rparam = static_cast<GenParam::RealParameter*>(getParameter(R_SOURCE));
    rparam->setMinValue(0.0);

    ParameterBase::GetVecFunc<Real> getFct = [&]()-> Real* { return m_boxMin.data(); };
	ParameterBase::SetVecFunc<Real> setFct = [&](Real* val)	
    {
	    m_boxMin = Vector3r(val[0], val[1], val[2]);
    };
	BOX_MIN = createVectorParameter("boxMin", "box min", 3u, getFct, setFct);
	setGroup(BOX_MIN, "coagualtion");
	setDescription(BOX_MIN, "Minimum point of box of which the rSource is not zero.");

	ParameterBase::GetVecFunc<Real> getFct2 = [&]()-> Real* { return m_boxMax.data(); };
    ParameterBase::SetVecFunc<Real> setFct2 = [&](Real* val)
    {
        m_boxMax = Vector3r(val[0], val[1], val[2]);
    };
	BOX_MAX = createVectorParameter("boxMax", "box max", 3u, getFct2, setFct2);
	setGroup(BOX_MAX, "coagualtion");
	setDescription(BOX_MAX, "Maximum point of box of which the rSource is not zero.");

    POINT_SRC_VAL = createNumericParameter("pointSrcVal", "pointSrcVal", &m_pointSrcVal);
	setDescription(POINT_SRC_VAL, "give a particle initial non-zero value to test the diffusion");
    setGroup(POINT_SRC_VAL, "coagualtion");
    rparam = static_cast<GenParam::RealParameter*>(getParameter(POINT_SRC_VAL));
    rparam->setMinValue(0.0);

    POINT_SRC_POS = createNumericParameter("pointSrcPos", "pointSrcPos", &m_pointSrcPos);
	setDescription(POINT_SRC_POS, "give a particle ID(uint) of which the initial value is not zero to test the diffusion");
    setGroup(POINT_SRC_POS, "coagualtion");

}

void Coagulation::step()
{
    Simulation* sim = Simulation::getCurrent();
    const unsigned int numParticles = m_model->numActiveParticles();
    const Real radius = sim->getValue<Real>(Simulation::PARTICLE_RADIUS);
    const Real diameter = static_cast<Real>(2.0) * radius;
    const Real diameter2 = diameter * diameter;
    const unsigned int fluidModelIndex = m_model->getPointSetIndex();
    const unsigned int nFluids = sim->numberOfFluidModels();
    const unsigned int nBoundaries = sim->numberOfBoundaryModels();
    const Real density0 = m_model->getDensity0();
    FluidModel* model = m_model;
    const Real h = sim->getSupportRadius();
    const Real dt = TimeManager::getCurrent()->getTimeStepSize();

    if(m_pointSrcPos!=-1)
        m_ccf[m_pointSrcPos] = m_pointSrcVal;

    // compute ccf
    #pragma omp parallel default(shared)
    {
        #pragma omp for schedule(static)  
        for (int i = 0; i < (int)numParticles; i++)
        {
            const Vector3r &xi = m_model->getPosition(i);
            Real &ccf_i = getCcf(i);
            Real density_i = m_model->getDensity(i);
            Real source = 0.0;
            Real ccf_sum = 0.0;
            Real ccf_old = ccf_i;
            
            if ((xi[0] > m_boxMin[0]) && (xi[1] > m_boxMin[1]) && (xi[2] > m_boxMin[2]) &&
				(xi[0] < m_boxMax[0]) && (xi[1] < m_boxMax[1]) && (xi[2] < m_boxMax[2]))
                source = m_rSource;
            for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i); j++)
            {
                const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, fluidModelIndex, i, j);
                const Vector3r &xj = model->getPosition(neighborIndex);
                const Real grawWNorm = sim->gradW(xi - xj).norm();

                Real density_j = m_model->getDensity(neighborIndex);
                Real ccf_j = getCcf(neighborIndex);
                ccf_sum += (m_diffusivity * m_model->getMass(neighborIndex) 
                / (density_j * density_i) * (ccf_j - ccf_i) * grawWNorm + source) * dt;
            }
            ccf_i = ccf_sum + ccf_old;
            model->setTemperature(i, ccf_i);
        }
    }

}

void Coagulation::reset()
{

}


void Coagulation::performNeighborhoodSearchSort()
{

    const unsigned int numPart = m_model->numActiveParticles();
    if (numPart == 0)
        return;

    Simulation* sim = Simulation::getCurrent();
    auto const& d = sim->getNeighborhoodSearch()->point_set(m_model->getPointSetIndex());
    d.sort_field(&m_ccf[0]);

}

void Coagulation::deferredInit()
{
    initValues();
}

void Coagulation::initValues()
{

}
