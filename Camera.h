#pragma once

#include <iostream>
#include <fstream>
#include <cmath>

#include <glm/glm.hpp>
#include "Ray.h"
#include "Transform.h"

using namespace std;

class Camera
{
	public:

	Camera() : eyepos(0.0f, 0.0f, 0.0f),
                forward(0.0f, 0.0f, 1.0f),
                up(0.0f, 1.0f, 0.0f),
                right(1.0f, 0.0f, 0.0f),
                vfov(DEG_TO_RAD * 60.0f)
    {
        aspect_ratio = 640.f / 480.f;
        calculateBottomLeft();
    }
	~Camera(){ }

    void SetResolution(uint32_t width, uint32_t height)
    {
        aspect_ratio = (float)width / (float)height;
        calculateBottomLeft();
    }

    void SetFOV(float theta) // horizontal field of view in degrees
    {
        vfov = DEG_TO_RAD * theta;
        calculateBottomLeft();
    }

    Ray generateRay(float u, float v)
    {
        return Ray(eyepos, bottomleft + u*horizontal_size + v*vertical_size);
    }

	public:
    float vfov; //vertical field of view
    float aspect_ratio;
    glm::vec3 eyepos;
	glm::vec3 forward;
	glm::vec3 up;
	glm::vec3 right;
    Transform transform;

    private:
    glm::vec3 bottomleft;
    glm::vec3 horizontal_size;
    glm::vec3 vertical_size;
    
    // This function shoud be called if eyepos, vfov, width, or height are changed
    void calculateBottomLeft()
    {
        // Assume the distance to the viewing plane is 1, and solve for width
        float viewport_height = 2.0f * tan(0.5f * vfov);
        float viewport_width = aspect_ratio * viewport_height;

        horizontal_size = glm::vec3(viewport_width, 0.f, 0.f);
        vertical_size = glm::vec3(0.f, viewport_height, 0.f);

        bottomleft = eyepos + forward - 0.5f * horizontal_size - 0.5f * vertical_size;

        std::cout << "Bottom left is " << bottomleft.x  << "," << bottomleft.y << "," << bottomleft.z << std::endl;
    }
};