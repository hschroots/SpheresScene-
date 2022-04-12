#pragma once

#include <glm/glm.hpp>

struct Light
{
	Light(){}
	~Light(){}
	glm::vec3 position;
	glm::vec3 color;
};

struct PointLight
{
    PointLight(){}
	~PointLight(){}
    glm::vec3 uniformDirection;
	glm::vec3 color;
};