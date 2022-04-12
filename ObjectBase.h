#pragma once

#include <glm/glm.hpp>
#include "Transform.h"

enum ShadingModel
{
    BlinnPhong = 0,
    ReflectOnly,
    RefractOnly,
    ReflectAndRefract,
    NUM_SHADING_MODELS
};

class ObjectBase
{
    public:
	ObjectBase(){};
	virtual ~ObjectBase() = 0;

    public:
    ShadingModel shadingType;
    glm::vec3 color;
	float IndexOfRefraction;
    Transform transform;
};
ObjectBase::~ObjectBase() {}