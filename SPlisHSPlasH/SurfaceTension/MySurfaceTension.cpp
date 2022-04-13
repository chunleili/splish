#include "MySurfaceTension.h"
#include "SPlisHSPlasH\Simulation.h"

using namespace SPH;


int MySurfaceTension::THRES_HIGH  = -1;
int MySurfaceTension::THRES_LOW = -1;
int MySurfaceTension::DIFFUSIVITY = -1;
int MySurfaceTension::R_SOURCE = -1;



//MySurfaceTension::MySurfaceTension(FluidModel *model) :
//	NonPressureForceBase(model),
//	m_youngsModulus(static_cast<Real>(100000.0)),
//	m_poissonRatio(static_cast<Real>(0.3))
//{
//	m_fixedBoxMin.setZero();
//	m_fixedBoxMax.setZero();
//}

MySurfaceTension::MySurfaceTension(FluidModel* model) :
    SurfaceTensionBase(model), 
    m_ccf(),
    m_thresHigh(static_cast<Real>(1.0)),
    m_thresLow(static_cast<Real>(0.1)),
    m_diffusivity(static_cast<Real>(0.1)),
    m_rSource(static_cast<Real>(0.1))
{


    m_ccf.resize(model->numParticles(), 0.0);

    model->addField({ "ccf field", FieldType::Scalar, [&](const unsigned int i) -> Real* { return &m_ccf[i]; } });

}



//MySurfaceTension::MySurfaceTension(FluidModel *model) :
//    SurfaceTensionBase(model), m_ccf()
//{
//    CF_threshold_low = 0.1;
//    CF_threshold_high = 0.5;
//    diffusivity = 0.1 ;
//    R_source = 0.1;
//
//    m_ccf.resize(model->numParticles(), 0.0);
//    
//    model->addField({ "ccf field", FieldType::Scalar, [&](const unsigned int i) -> Real* { return &m_ccf[i]; } });
//
//}

MySurfaceTension::~MySurfaceTension(void)
{
    m_model->removeFieldByName("ccf field");
    m_ccf.clear();
}

void MySurfaceTension::initParameters()
{
    SurfaceTensionBase::initParameters();
}

void MySurfaceTension::step()
{
    Simulation* sim = Simulation::getCurrent();
    const unsigned int numParticles = m_model->numActiveParticles();
    const Real k = m_surfaceTension;
    const Real kb = m_surfaceTensionBoundary;
    const Real radius = sim->getValue<Real>(Simulation::PARTICLE_RADIUS);
    const Real diameter = static_cast<Real>(2.0) * radius;
    const Real diameter2 = diameter * diameter;
    const unsigned int fluidModelIndex = m_model->getPointSetIndex();
    const unsigned int nFluids = sim->numberOfFluidModels();
    const unsigned int nBoundaries = sim->numberOfBoundaryModels();
    const Real density0 = m_model->getDensity0();
    FluidModel* model = m_model;
    const Real h = sim->getSupportRadius();

    // Compute forces
//#pragma omp parallel default(shared)
    {
//#pragma omp for schedule(static)  
        for (int i = 0; i < (int)numParticles; i++)
        {
            const Vector3r &xi = m_model->getPosition(i);
            Real &ccf_i = getCcf(i);

            Real density_i = m_model->getDensity(i);
            //////////////////////////////////////////////////////////////////////////
            // Fluid
            //////////////////////////////////////////////////////////////////////////
            forall_fluid_neighbors_in_same_phase(
                Real density_j = m_model->getDensity(neighborIndex);
                Real ccf_j = getCcf(neighborIndex);
                ccf_i += m_diffusivity * m_model->getMass(neighborIndex) / (density_j * density_i) * (ccf_i - ccf_j) * sim->gradW(xi - xj).norm();
            );
        }
    }
}

void MySurfaceTension::reset()
{

}


void MySurfaceTension::performNeighborhoodSearchSort()
{

    const unsigned int numPart = m_model->numActiveParticles();
    if (numPart == 0)
        return;

    Simulation* sim = Simulation::getCurrent();
    auto const& d = sim->getNeighborhoodSearch()->point_set(m_model->getPointSetIndex());
    d.sort_field(&m_ccf[0]);

}

void MySurfaceTension::deferredInit()
{
    initValues();
}

void MySurfaceTension::initValues()
{
    //TODO: tune parameter in GUI (should be in the base class)

 //   THRES_LOW = createNumericParameter("THRES_LOW", "THRES_LOW", &m_thresLow);
	//setGroup(THRES_LOW, "SurfaceTension");
	//setDescription(THRES_LOW, "lower threshold of coagualtion factor");
	//RealParameter* rparam = static_cast<RealParameter*>(getParameter(THRES_LOW));
	//rparam->setMinValue(0.0);

 //   THRES_HIGH = createNumericParameter("THRES_HIGH", "THRES_HIGH", &m_thresHigh);
	//setGroup(THRES_HIGH, "SurfaceTension");
	//setDescription(THRES_HIGH, "higher threshold of coagualtion factor");
	//RealParameter* rparam = static_cast<RealParameter*>(getParameter(THRES_HIGH));
	//rparam->setMinValue(0.0);

 //   DIFFUSIVITY = createNumericParameter("DIFFUSIVITY", "DIFFUSIVITY", &m_diffusivity);
	//setGroup(DIFFUSIVITY, "SurfaceTension");
	//setDescription(DIFFUSIVITY, "higher threshold of coagualtion factor");
	//RealParameter* rparam = static_cast<RealParameter*>(getParameter(DIFFUSIVITY));
	//rparam->setMinValue(0.0);

 //   R_SOURCE = createNumericParameter("R_SOURCE", "R_SOURCE", &m_rSource);
	//setGroup(R_SOURCE, "SurfaceTension");
	//setDescription(R_SOURCE, "higher threshold of coagualtion factor");
	//RealParameter* rparam = static_cast<RealParameter*>(getParameter(R_SOURCE));
	//rparam->setMinValue(0.0);
}