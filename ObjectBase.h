#pragma once

#include <glm/glm.hpp>
#include "Transform.h"

enum ShadingModel
{
    Opaque = 0,
    Transparent,
    NUM_SHADING_MODELS
};

class ObjectBase
{
    public:
	ObjectBase(){};
	virtual ~ObjectBase() = 0;
    virtual glm::vec3 ColorAt(const glm::vec3 position) = 0;

    public:
    ShadingModel shadingType;
    glm::vec3 color;
    float reflectivity;
    float transmittance;
    float IndexOfRefraction;
    Transform transform;
};
ObjectBase::~ObjectBase() {}