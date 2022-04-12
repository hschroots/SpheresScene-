#pragma once

#include <numeric>
#include <glm/glm.hpp>

struct Ray
{
    Ray() : origin(glm::vec3(0.f, 0.f, 0.f)),
            direction(glm::vec3(0.f,0.f,1.f)),
            tmin(0.f),
            tmax(std::numeric_limits<float>::max())
    {
        //empty
    }
    Ray(glm::vec3 nOrigin, glm::vec3 nDirection) : 
        origin(nOrigin),
        direction(nDirection),
        tmin(0.f),
        tmax(std::numeric_limits<float>::max())
    {
        //empty
    }
    ~Ray() {}

    inline glm::vec3 at(float& t) const
    {
        return origin + t * direction;
    }
    glm::vec3 origin;
    glm::vec3 direction;
    float tmin;
    float tmax;
};