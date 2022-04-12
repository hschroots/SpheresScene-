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
};