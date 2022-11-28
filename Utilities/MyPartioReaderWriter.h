#pragma once

#include "SPlisHSPlasH/Common.h"
#include <vector>

namespace Utilities
{

	class MyPartioReaderWriter
	{
	public:
		static bool readParticles(const std::string &fileName, std::vector<Vector3r> &pos);
	};

}

