#include "MySurfaceTension.h"

namespace SPH{

MySurfaceTension::MySurfaceTension(FluidModel *model) :
	SurfaceTensionBase(model)
{

}

MySurfaceTension::~MySurfaceTension(void)
{

}

void MySurfaceTension::initParameters()
{
	SurfaceTensionBase::initParameters();
}

void MySurfaceTension::step()
{

}

void MySurfaceTension::reset()
{

}

}