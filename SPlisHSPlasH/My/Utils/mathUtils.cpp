#include <vector>
#include "mathUtils.h"

namespace SPH
{

    //对一个场求和
    float sumField(std::vector<float>& field, unsigned int numParticles)
    {
        float sum = 0.0f;
        for (unsigned int i = 0; i < numParticles; i++)
        {
            sum += field[i];
        }
        return sum;
    }

    float averageField(std::vector<float>& field, unsigned int numParticles)
    {
        return sumField(field, numParticles) / numParticles;
    }

    float maxField(std::vector<float>& field, unsigned int numParticles)
    {
        float max = field[0];
        for (unsigned int i = 1; i < numParticles; i++)
        {
            if (field[i] > max)
                max = field[i];
        }
        return max;
    }

    float minField(std::vector<float>& field, unsigned int numParticles)
    {
        float min = field[0];
        for (unsigned int i = 1; i < numParticles; i++)
        {
            if (field[i] < min)
                min = field[i];
        }
        return min;
    }
    
} // namespace SPH
