#pragma once

#include "SPlisHSPlasH/Common.h"
#include <vector>

namespace Utilities
{

	class MyPartioReader
	{
	public:
		static bool readParticles(const std::string &fileName, std::vector<Vector3r> &pos);
		static bool  readParticlesUv(const std::string &fileName, std::vector<Vector3r> &uv);
	};

}

