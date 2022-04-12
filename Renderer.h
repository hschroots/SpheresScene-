#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <exception>

#include "ObjectBase.h"
#include "Sphere.h"
#include "Plane.h"
#include "Camera.h"
#include "Image.h"
#include "Ray.h"

using namespace std;

class Renderer
{
	public:

	Renderer(){ }
	~Renderer(){ }

	void setFramebufferDimensions(uint32_t wWidth = 640, uint32_t wHeight = 480)
	{
		width = wWidth;
		height = wHeight;
		framebuffer = Image(width, height);
		camera.SetResolution(width, height);
	}
	void init()
	{
		setFramebufferDimensions();
		camera.SetFOV(60.0f); //horizontal field of view in degrees
	}

	string printVec(glm::vec3 v) const
    {
		stringstream ss;
        ss << "(" << v.x << "," << v.y << "," << v.z << ")";
		return ss.str();
    }

	bool intersectSphere(Ray& ray, Sphere obj, HitRecord& hitRec, bool debug)
	{
		//float t_imageplane = (float)glm::
		glm::vec3 uray = glm::normalize(ray.direction);
		//std::cout << "uRay: " << uray.x << " , " << uray.y << ", " << uray.z << std::endl;

		float eps = 0.01f;
		float rr = obj.Radius() * obj.Radius();
		glm::vec3 diff = obj.Center() - ray.origin;
		float t0 = glm::dot(diff, uray);
		float dSquared = glm::dot(diff, diff) - t0 * t0;

		if(debug)
		{
			std::cout << "uray" << printVec(uray) << std::endl;
			std::cout << "obj.center" << printVec(obj.Center()) << std::endl;
			std::cout << "obj.radiud^2 = " << rr << std::endl;
			std::cout << "Vec from eye to obj.center(diff) " << printVec(diff) << std::endl;
			std::cout << "t0 = dot(diff, uray) = " << t0 << std::endl;
			std::cout << "glm::dot(diff, diff) = " << glm::dot(diff, diff) << std::endl;
			std::cout << "dSquared = glm::dot(diff, diff) - t0 * t0 = " << dSquared << std::endl;
		}

		if( dSquared > rr)
		{
			return false;
		}

		float t1 = sqrt( rr - dSquared );
		float intersectionDistance = t0 > t1 + eps ? t0 - t1 : t0 + t1;

		if(debug)
		{
			std::cout << "t1 = sqrt( rr - dSquared ) = " << t1 << std::endl;
			std::cout << "intersection distance =" << intersectionDistance << std::endl;
		}
		if (intersectionDistance > eps) // we have a hit!
		{
			hitRec.hitPoint = ray.at(intersectionDistance);
			hitRec.hitNormal = (hitRec.hitPoint - obj.Center()); // / obj.raidus;
			hitRec.hitNormal = glm::normalize(hitRec.hitNormal);

			if(debug)
			{
				std::cout << "pHit = " << printVec(hitRec.hitPoint) << std::endl;
				std::cout << "nHit = " << printVec(hitRec.hitNormal) << std::endl;
			}
			return true;
		}
		return false;
	}
	bool intersectPlane(glm::vec3 rayOrigin, glm::vec3 ray, Plane obj, HitRecord& rec, bool debug)
	{
		glm::vec3 uray = glm::normalize(ray);
		glm::vec3 n = obj.planeNormal;
		float eps = 0.000001;
		// assuming vectors are all normalized
		float denom = glm::dot(n, uray);
		if (std::abs(denom) > eps) { 
			glm::vec3 p0l0 = obj.planePosition - rayOrigin; 
			//p0l0 = glm::normalize(p0l0);
			float t = glm::dot(p0l0, n) / denom;
			if(t >= 0)
			{
				rec.numHits = 1;
				rec.t = t;
				rec.hitPoint = rayOrigin + (t * uray);
				rec.hitNormal = glm::normalize(n);
				return true;
			}
		} 
		return false;
	}

	bool intersect(Ray& ray, ObjectBase* obj, HitRecord& rec, bool debug)
	{
		// int dataType
		// 0 - Sphere
		// 1 - Plane

		int dataType = -1;
		
		if (dynamic_cast<Sphere*>(obj) != nullptr)
		{
			//std::cout << "Not a sphere" << std::endl;
			dataType = 0;
		}
		else if (dynamic_cast<Plane*>(obj) != nullptr)
		{
			//std::cout << "Not a plane" << std::endl;
			dataType = 1;
		}

		switch(dataType)
		{
			case 0:
					{
					Sphere *s = dynamic_cast<Sphere*>(obj);
					if(s->intersect(ray, rec, debug))
					//if(intersectSphere(ray, *s, hitRec, debug))
					{
						return true;
					}
					return false;
					}
			break;
			case 1:
					{
					Plane *p = dynamic_cast<Plane*>(obj);
					return intersectPlane(ray.origin, ray.direction, *p, rec, debug);
					}
			break;
			default:
				return false;
		}
	}

	void computeRay(float u, float v, Ray& ray)
    {
        ray = camera.generateRay(u, v);
    }

	public:
    Camera camera;
	Image framebuffer;
	uint32_t width;
	uint32_t height;
};