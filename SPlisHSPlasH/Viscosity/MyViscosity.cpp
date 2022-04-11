#include "SPlisHSPlasH\Simulation.h"
#include "MyViscosity.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "SPlisHSPlasH/BoundaryModel_Koschier2017.h"
#include "SPlisHSPlasH/BoundaryModel_Bender2019.h"

//#include "extern\cuNSearch\src\Ext_NeighborhoodSearch\include\PointSet.h"


using namespace SPH;
using namespace GenParam;

MyViscosity::MyViscosity(FluidModel *model) :
	ViscosityBase(model), m_vDiff()
{
	m_vDiff.resize(model->numParticles(), Vector3r::Zero());
	
	model->addField({ "MyAddedField", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_vDiff[i][0]; }, true });
}

MyViscosity::~MyViscosity(void)
{
	m_model->removeFieldByName("MyAddedField");
	m_vDiff.clear();
}

void MyViscosity::initParameters()
{
	ViscosityBase::initParameters();


}

void MyViscosity::step()
{

}

void MyViscosity::reset()
{

}

void MyViscosity::performNeighborhoodSearchSort()
{
	//Simulation* sim = Simulation::getCurrent();
	//auto const& d = sim->getNeighborhoodSearch()->point_set(m_model->getPointSetIndex());
	//d.sort_field(&m_myParticleViscosityData[0]);


	const unsigned int numPart = m_model->numActiveParticles();
	if (numPart == 0)
		return;

	Simulation* sim = Simulation::getCurrent();
	auto const& d = sim->getNeighborhoodSearch()->point_set(m_model->getPointSetIndex());
	d.sort_field(&m_vDiff[0]);

}

void MyViscosity::deferredInit()
{
	initValues();
}

void MyViscosity::initValues()
{
}