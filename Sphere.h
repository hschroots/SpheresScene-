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
	inline glm::vec3 NormalAt(glm::vec3& worldPoint) const
	{
		//convert worldpoint to object point
		glm::vec3 objectPoint = transform.WorldToObject() * glm::vec4(worldPoint, 1.f);
		glm::vec3 objectNormal = (objectPoint - center) / radius;
		// Convert object space nromal to world space by the inverse transpose of the model matrix
		glm::mat4x4 invObject = glm::inverse(transform.ObjectToWorld());
		glm::mat4x4 invTrans = glm::transpose(invObject);
		glm::vec3 worldNormal = invTrans * glm::vec4(objectNormal, 0.f);
		return glm::normalize(worldNormal);
	}

	glm::vec3 ColorAt(const glm::vec3 position) override
	{
		return color;
	}

	bool intersect(Ray& ray, HitRecord& hitRec, string type, bool debug)
	{
		// Assume ray is normalized
		// Convert ray to object space
		Ray objectRay;
		objectRay.origin	 = transform.WorldToObject() * glm::vec4(ray.origin   , 1.f);
		objectRay.direction  = transform.WorldToObject() * glm::vec4(ray.direction, 0.f);

		glm::vec3 sphere_to_ray = objectRay.origin - center;
		float directionMagnitude = glm::length(objectRay.direction);
		float ocMagnitude = glm::length(sphere_to_ray);
		// Do all calculations in object space
		float a = glm::dot(objectRay.direction, objectRay.direction); // dot product of self is same as magnitude sqaured
		float half_b = glm::dot(sphere_to_ray, objectRay.direction);
		float c = glm::dot(sphere_to_ray, sphere_to_ray) - radius * radius;
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
			std::cout << type + " " << "RayOrigin" << printVec(ray.origin) << std::endl;
			std::cout << type + " " << "ObjectRayOrigin" << printVec(objectRay.origin) << std::endl;
			std::cout << type + " " << "rayDir" << printVec(ray.direction) << std::endl;
			std::cout << type + " " << "objectRayDir" << printVec(objectRay.direction) << std::endl;
			std::cout << type + " " << "sphere center" << printVec(center) << std::endl;
			std::cout << type + " " << "sphere radius = " << radius << std::endl;
			std::cout << std::flush;

			/*
			std::cout << "Vec from eye to obj.center(diff) " << printVec(diff, std::cout) << std::endl;
			std::cout << "t0 = dot(diff, uray) = " << t0 << std::endl;
			std::cout << "glm::dot(diff, diff) = " << glm::dot(diff, diff) << std::endl;
			*/
			std::cout << type + " " << "dSquared =  " << dSquared << std::endl;
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
			//std::cout << type + " " << "t0 outside valid range" << std::endl;
			intersectionDistance = t1;
			if (t1 < objectRay.tmin || t1 > objectRay.tmax)
			{
				//std::cout << type + " " << "t1 outside valid range" << std::endl;
				// both values are outisde of the value range, no hit
				return false;
			}
		}

		if (abs(t1 - t0) < eps)
			hitRec.numHits = 1;
		else
			hitRec.numHits = 2;

		if(debug)
		{
			std::cout << type + " " << "t0 = " << t0  << " t1 = "  << t1 << std::endl;
			std::cout << type + " " << "intersection distance = " << intersectionDistance << std::endl;
		}
		if (intersectionDistance > eps) // we have a hit!
		{
			hitRec.t = intersectionDistance;
			glm::vec3 objectHitPoint = objectRay.at(intersectionDistance);
			hitRec.hitPoint = transform.ObjectToWorld() * glm::vec4(objectHitPoint, 1.f);
			glm::vec3 objectNormal = (objectHitPoint - center) / radius;
			if(debug)
			{
				std::cout << type + " " << "objectHitPoint = " << printVec(objectHitPoint) << std::endl;
				std::cout << type + " " << "hitPoint = " << printVec(hitRec.hitPoint) << std::endl;
				std::cout << type + " " << "distance b/t the two = " << glm::distance(objectHitPoint,hitRec.hitPoint) << std::endl;
				std::cout << type + " " << "objectNormal = " << printVec(objectNormal) << std::endl;
				if(glm::length(objectNormal) < 0.99999f)
					std::cout << "objectNormal has length " << glm::length(objectNormal) << std::endl;
			}

	
			//If the ray direction and the normal are negative then we're a front facing ray hitting the sphere surface
			if (glm::dot(objectRay.direction, objectNormal) <= 0.f)
				hitRec.front_face = true;
			else // our ray is going in the same direction as the normal, so we must be inside the sphere
			{
				hitRec.front_face = false;
				//objectNormal = -objectNormal;
			}
			
			// Tranform normals by the inverse transpose
			glm::mat4x4 invObject = glm::inverse(transform.ObjectToWorld());
			glm::mat4x4 invTrans = glm::transpose(invObject);
			glm::vec3 tmpNorm = invTrans * glm::vec4(objectNormal, 0.f);
			hitRec.hitNormal = glm::normalize(tmpNorm);

			if(debug)
			{
				std::cout << type + " " << "Is front face = " << (hitRec.front_face ? "true" : "false") << std::endl;
				std::cout << type + " " << "hitNormal = " << printVec(hitRec.hitNormal) << std::endl;
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