#pragma once

#include <vector>
#include "SPlisHSPlasH/FluidModel.h"
#include "SPlisHSPlasH/Simulation.h"

//获取哪些粒子是表面粒子
namespace SPH
{
    
void findSurfaceParticles(FluidModel *model, std::vector<int> &id)
{
    Simulation* sim = Simulation::getCurrent();
    const int numParticles = model->numActiveParticles();
    const unsigned int fluidModelIndex =  model->getPointSetIndex();
    for (int i = 0; i < numParticles; i++)
    {
        int numNeighbors = sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i);
        // 如果粒子的邻居数小于10，那么就是表面粒子
        if (numNeighbors < 10)
        {
            id.push_back(i);
        }
    }
}

} // namespace SPH