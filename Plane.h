#pragma once

#include <glm/glm.hpp>
#include "ObjectBase.h"

class Plane : public ObjectBase
{
	public:
	Plane(){}
	~Plane(){}
    glm::vec3 planePosition;
	glm::vec3 planeNormal;
	bool useColorPattern;

	glm::vec3 ColorAt(glm::vec3 position) override
	{
		if(useColorPattern)
		{
			glm::vec3 dark = 0.5f * color;
			glm::vec3 light = color;

			float xCoord = floor(position.x);
			bool XisEven = (int)xCoord % 2 == 0;

			float zCoord = floor(position.z);
			bool ZisEven = (int)zCoord % 2 == 0;

			if( (XisEven && ZisEven) || (!XisEven && !ZisEven))
				return light;
			return dark;
		}
		else
			return color;
	}
};