#pragma once

#include <iostream>
#include <glm/glm.hpp>

#include "Ray.h"
#include "HitRecord.h"
#include "ObjectBase.h"

using namespace std;

class Sphere : public ObjectBase
{
	public:
	Sphere() : center(glm::vec3(0.f)),
				radius(1.f)
	{
		//empty
	}
	~Sphere(){}

	//GetMethods
	inline glm::vec3 Center() const { return transform.ObjectToWorld() * glm::vec4(center, 1.f); }
	inline float Radius() const { return radius; }
	inline glm::vec3 NormalAt(glm::vec3& point) const
	{
		//convert worldpoint to object point
		glm::vec3 objectPoint = transform.WorldToObject() * glm::vec4(point, 1.f);
		glm::vec3 objectNormal = objectPoint - center;
		return glm::transpose(transform.WorldToObject()) * glm::vec4(objectNormal, 0.f);
	}

	bool intersect(Ray& ray, HitRecord& hitRec, bool debug)
	{
		// Assume ray is normalized
		// Convert ray to object space
		Ray objectRay = ray;
		objectRay.origin = transform.WorldToObject() * glm::vec4(ray.origin, 1.f);
		objectRay.direction    = transform.WorldToObject() * glm::vec4(ray.direction, 0.f);

		glm::vec3 sphere_to_ray = objectRay.origin - center;
		float directionMagnitude = glm::length(objectRay.direction);
		float ocMagnitude = glm::length(sphere_to_ray);
		// Do all calculations in object space
		float a = directionMagnitude * directionMagnitude; // dot product of self is same as magnitude sqaured
		float half_b = glm::dot(objectRay.direction, sphere_to_ray);
		float c = ocMagnitude * ocMagnitude - radius * radius;
		float dSquared = half_b * half_b - a * c;
		
		float eps = 0.00000001f;
		/*
		float rr = radius * radius;
		glm::vec3 diff = center - objectRayOrigin;
		float t0 = glm::dot(diff, objectRayDir);
		float dSquared = glm::dot(diff, diff) - t0 * t0;
		*/
		if(debug)
		{
			std::cout << "RayOrigin" << printVec(objectRay.origin) << std::endl;
			std::cout << "ObjectRayOrigin" << printVec(objectRay.origin) << std::endl;
			std::cout << "ray" << printVec(ray.direction) << std::endl;
			std::cout << "objectRayDir" << printVec(objectRay.direction) << std::endl;
			std::cout << "sphere center" << printVec(center) << std::endl;
			std::cout << "sphere radius = " << radius * radius << std::endl;

			/*
			std::cout << "Vec from eye to obj.center(diff) " << printVec(diff, std::cout) << std::endl;
			std::cout << "t0 = dot(diff, uray) = " << t0 << std::endl;
			std::cout << "glm::dot(diff, diff) = " << glm::dot(diff, diff) << std::endl;
			*/
			std::cout << "dSquared =  " << dSquared << std::endl;
		}
		
		if( dSquared < 0.f)
		{
			return false;
		}

		float d = sqrt(dSquared);
		float t0 = ( -half_b - d ) / a;
		float t1 = ( -half_b + d ) / a;
		float intersectionDistance = t0;
		if (t0 < objectRay.tmin || t0 > objectRay.tmax)
		{
			intersectionDistance = t1;
			if (t1 < objectRay.tmin || t1 > objectRay.tmax)
				// both values are outisde of the value range, no hit
				return false;
		}

		if (abs(t1 - t0) < eps)
			hitRec.numHits = 1;
		else
			hitRec.numHits = 2;

		if(debug)
		{
			std::cout << "t0 = " << t0  << " t1 = "  << t1 << std::endl;
			std::cout << "intersection distance =" << intersectionDistance << std::endl;
		}
		if (intersectionDistance > eps) // we have a hit!
		{
			hitRec.t = intersectionDistance;
			glm::vec3 objectHitPoint = objectRay.at(intersectionDistance);
			hitRec.hitPoint = transform.ObjectToWorld() * glm::vec4(objectHitPoint, 1.f);
			glm::vec3 objectNormal = objectHitPoint - center;

			//If the ray direction and the normal are negative then we're a front facing ray hitting the sphere surface
			if (glm::dot(objectRay.direction, objectNormal) <= 0.f)
				hitRec.front_face = true;
			else // our ray is going in the same direction as the normal, so we must be inside the sphere
			{
				hitRec.front_face = false;
				objectNormal = -objectNormal;
			}
			// Tranform normals by the inverse transpose
			hitRec.hitNormal = glm::normalize(glm::transpose(transform.WorldToObject()) * glm::vec4(objectNormal, 0.f));

			if(debug)
			{
				std::cout << "hitPoint = " << printVec(hitRec.hitPoint) << std::endl;
				std::cout << "hitNormal = " << printVec(hitRec.hitNormal) << std::endl;
			}
			return true;
		}
		return false;
	}
	
	private:
	string printVec(glm::vec3 v) const
    {
		stringstream ss;
		ss << "(" << v.x << "," << v.y << "," << v.z << ")";
		return ss.str();
    }

	private:
	glm::vec3 center;
	float radius;
};