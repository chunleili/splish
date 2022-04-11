#include "MyVorticity.h"

using namespace SPH;
using namespace GenParam;

int MyVorticity::VORTICITY_COEFFICIENT = -1;

MyVorticity::MyVorticity(FluidModel *model) :
	VorticityBase(model)
{
	m_vorticityCoeff = static_cast<Real>(0.01);
}

MyVorticity::~MyVorticity(void)
{
}

void MyVorticity::initParameters()
{
	VorticityBase::initParameters();

	VORTICITY_COEFFICIENT = createNumericParameter("vorticity", "Vorticity transfer coefficient", &m_vorticityCoeff);
	setGroup(VORTICITY_COEFFICIENT, "Vorticity");
	setDescription(VORTICITY_COEFFICIENT, "Coefficient for the vorticity force computation");
	RealParameter* rparam = static_cast<RealParameter*>(getParameter(VORTICITY_COEFFICIENT));
	rparam->setMinValue(0.0);
}

void MyVorticity::step()
{

}

void MyVorticity::reset()
{
	
}

void MyVorticity::performNeighborhoodSearchSort()
{

}
