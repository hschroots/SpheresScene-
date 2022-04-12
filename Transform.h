#pragma once

#include <glm/glm.hpp>
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtc/matrix_inverse.hpp"

#define PI 3.14159265358979323846f
#define DEG_TO_RAD PI/180.f

class Transform
{
    public:
    Transform() :
        position(glm::vec3(0.f)),
        scale(glm::vec3(1.f)),
        rotation(glm::vec3(0.f)),
        objectToWorld(glm::mat4x4(1.f)),
        worldToObject(glm::mat4x4(1.f))
    {
        //empty
        hasChanged = false;
    }
    ~Transform(){}

    // All the "Get" methods
    inline glm::mat4x4 ObjectToWorld()const { return objectToWorld; }
    inline glm::mat4x4 WorldToObject()const{ return worldToObject; }
    inline glm::vec3 Position()const { return position; }
    inline glm::vec3 Scale()const { return scale; }
    //inline glm::vec3 EulerAngles()const { return rotation; }

    void Translate(glm::vec3& translate)
    {
        hasChanged = true;
        objectToWorld  = glm::translate(objectToWorld, translate);
        position = objectToWorld * glm::vec4(position,1.f);
        worldToObject = glm::inverse(objectToWorld);
        hasChanged = false;
    }
    void Scale(glm::vec3& scale)
    {
        hasChanged = true;
        objectToWorld = glm::scale(objectToWorld, scale);
        float xscale = objectToWorld[0][0];
        float yscale = objectToWorld[1][1];
        float zscale = objectToWorld[2][2];
        scale = glm::vec3(xscale, yscale, zscale);
        worldToObject = glm::inverse(objectToWorld);
        hasChanged = false;
    }

    // Roation in degrees
    void Rotate(float& degrees, glm::vec3& axis)
    {
        hasChanged = true;
        float radians = DEG_TO_RAD * degrees;
        objectToWorld = glm::rotate(objectToWorld, radians, axis);
        worldToObject = glm::inverse(objectToWorld);
        hasChanged = false;
    }

    private:
    bool hasChanged;
    glm::vec3 position;
    glm::vec3 scale;
    glm::vec3 rotation;
    glm::mat4x4 objectToWorld;
    glm::mat4x4 worldToObject;
};