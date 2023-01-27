#pragma once
#include <vector>

namespace SPH
{
    float sumField(std::vector<float>& field, unsigned int numParticles);
    float averageField(std::vector<float>& field, unsigned int numParticles);
    float maxField(std::vector<float>& field, unsigned int numParticles);
    float minField(std::vector<float>& field, unsigned int numParticles);
    auto minMaxField(std::vector<float>& field, unsigned int numParticles);
    
} // namespace SPH
