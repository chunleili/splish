#include "Partio.h"
using namespace Partio;
int main()
{
    //读入文件并构造ParticlesDataMutable对象data
    Partio::ParticlesDataMutable* data=Partio::read("E:/Dev/SPlisHSPlasH/extern/my_partio/data/cruiser.bgeo.sc");

    std::cout<<data->numParticles()<<"\n";

    Partio::write("E:/codes/try/partio/my_partio/data/cruiser_out.bgeo.sc",*data);

    // 释放
    data->release();
}