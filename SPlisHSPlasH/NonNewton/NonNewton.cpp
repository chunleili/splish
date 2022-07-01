#include "NonNewton.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/TimeManager.h"
using namespace SPH;

int NonNewton::NON_NEWTON = -1;

NonNewton::NonNewton(FluidModel *model) :
NonPressureForceBase(model),
m_viscosity(0.0)
{
	std::cout<<"constructor\n";
}

NonNewton::~NonNewton(void)
{
}

void NonNewton::init()
{
	std::cout<<"init\n";
	initParameters();
}

void NonNewton::step()
{
	static int steps{0};
	std::cout<<"step: "<<steps<<"\t";
	std::cout<<"m_viscosity: "<<m_viscosity<<"\n";
	steps++;
	m_viscosity++;
}

void NonNewton::reset()
{
	std::cout<<"reset\n";
	m_viscosity=0.0;
}

void NonNewton::initParameters()
{
	NON_NEWTON = createNumericParameter("nonNewtonViscosity", "nonNewton viscosity", &m_viscosity);
	setGroup(NON_NEWTON, "NonNewton");
	setDescription(NON_NEWTON, "NonNewton Params");
}