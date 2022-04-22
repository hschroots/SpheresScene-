#include "SpheresScene.h"

#include <cassert>

using namespace std;

string printVec(glm::vec3 v)
{
	stringstream ss;
	ss << "(" << v.x << "," << v.y << "," << v.z << ")";
	return ss.str();
}

void TestRaySphereIntersection()
{
	// There are multiple scenarios that we need to consider
	// Teh most common ones:
	// Family 1: The ray's origin is outside the sphere
	// case 1: the ray completely missed the sphere
	// case 2: the ray intersects the sphere trangentally
	// case 3: the ray intersects the sphere 
	// 
	// Family 2: the ray's origin is inside the sphere
	// The ray is guaranteed to intersect with the sphere
	// at one point in front of the ray and on behind it
	// the surface normal will be wrong

	Sphere s; // create a unit sphere at the origin

	// Create a ray whose origin is -5 on the z axis and +2 on the y-axis
	// This ray is parallel to the x-z plane.
	// Becuase the sphere is a radius of 1, the max y of any point in the sphere will be 1
	// So this ray should completely miss sphere
	Ray ray;
	ray.origin = glm::vec3(0.f, 2.0f, -5.f);
	ray.direction = glm::vec3(0.f, 0.f, 1.f);

	std::cout << " --- BEGIN TEST ---" << std::endl;
	std::cout << " --- Test Ray missing intersection with sphere ---" << std::endl;
	{
		HitRecord rec;
		bool hit = s.intersect(ray, rec, "P", false);
		assert(hit != true);
	}

	// Now create a ray this is exactly +1  on the y axis
	// This ray should intersect the sphere tangentailly
	ray.origin = glm::vec3(0.f, 1.f, -5.0f);

	std::cout << " --- Test Ray intersecting tangentially with spehere (1 root) ---" << std::endl;
	{
		HitRecord rec;
		bool hit = s.intersect(ray, rec, "P",false);
		assert(hit == true);
		assert(rec.numHits == 1);
		assert(rec.front_face == true);

		// Check that the normal has been normalized
		assert(glm::length(rec.hitNormal) >= 0.9999f);
		//Check that the ray direction and the sphere normal are negative
		//assert(glm::dot(ray.direction, rec.hitNormal) <= 0.f);
	}

	// Now create a ray that will interserct the sphere twice
	ray.origin = glm::vec3(0.f, 0.0f, -5.0f);
	std::cout << " --- Test Ray full intersection with sphere ---" << std::endl;
	{
		HitRecord rec;
		bool hit = s.intersect(ray, rec, "P",true);
		assert(hit == true);
		assert(rec.numHits == 2);
		assert(rec.front_face == true);

		// Check that the normal has been normalized
		assert(glm::length(rec.hitNormal) >= 0.9999f);
		//Check that the ray direction and the sphere normal are negative
		//assert(glm::dot(ray.direction, rec.hitNormal) <= 0.f);
	}

	// Now create a ray whos' origin is inside the sphere
	ray.origin = glm::vec3(0.f, 0.f, 0.f);
	std::cout << " --- Test Ray orgin inside sphere ---" << std::endl;
	{
		HitRecord rec;
		bool hit = s.intersect(ray, rec, "P",true);
		assert(hit == true);
		assert(rec.numHits == 2);
		assert(rec.front_face == false);

		// Check that the normal has been normalized
		assert(glm::length(rec.hitNormal) >= 0.9999f);
		//Check that the ray direction and the sphere normal are negative
		//assert(glm::dot(ray.direction, rec.hitNormal) <= 0.f);
	}
	std::cout << " --- Test Ray orgin inside sphere ---" << std::endl;
	std::cout << " --- END TEST ---" << std::endl;
}

void TestReflectDirection()
{
	// Dealing with floats so check against epsilon, not 0
	float eps2 = std::numeric_limits<float>::epsilon() * 2.0f;
	std::cout << " --- BEGIN REFLECT TEST ---" << std::endl;
	glm::vec3 ray = glm::vec3(0.f, -1.f, -1.f);
	ray = glm::normalize(ray); // this should be a ray downward at 45 degrees

	glm::vec3 normal = glm::vec3(0.f, 1.f, 0.f);
	
	std::cout << " --- Reflecting Ray ---" << std::endl;
	glm::vec3 reflectedRay = glm::reflect(ray,normal);

	//Check that the ray has lenght 1
	std::cout << "Ray length after reflections is " << glm::length(reflectedRay) << std::endl << std::flush;
	assert(glm::length(reflectedRay) > 0.99999f );

	//Check that the ray has been reflected
	glm::vec3 expectedReflectedRay = glm::vec3(0.f, 1.f, -1.f);
	expectedReflectedRay = glm::normalize(expectedReflectedRay);

	for(int i = 0; i < expectedReflectedRay.length(); i++)
	{
		assert( abs(reflectedRay.x - expectedReflectedRay.x) < eps2);
		assert( abs(reflectedRay.y - expectedReflectedRay.y) < eps2);
		assert( abs(reflectedRay.z - expectedReflectedRay.z) < eps2);
	}
	std::cout << " --- Ray Reflected as expected ---" << std::endl;
	ray = glm::vec3(-1.f, -1.f, -1.f); // Ray should be at 45 degrees in all 3 planes
	ray = -1.f * glm::normalize(ray);

	std::cout << " --- Reflecting negated Ray ---" << std::endl;
	reflectedRay = glm::reflect(ray,normal);

	for(int i = 0; i < expectedReflectedRay.length(); i++)
	{
		assert( false == (abs(reflectedRay.x - expectedReflectedRay.x) < eps2));
		assert( false == (abs(reflectedRay.y - expectedReflectedRay.y) < eps2));
		assert( false == (abs(reflectedRay.z - expectedReflectedRay.z) < eps2));
	}
	std::cout << " --- END REFLECT TEST ---" << std::endl;
}
void TestFresnelNumber()
{
	float eps = std::numeric_limits<float>::epsilon();
	float eps2 = eps * 2.0f;
	std::cout << " --- BEGIN FRESNEL TEST ---" << std::endl;

	// Create a sphere with refractive index near that of glass
	Sphere s;
	s.IndexOfRefraction = 1.5;

	// Create a ray that will intersect with the sphere at a 45 degree angle
	Ray incidentRay;
	incidentRay.origin = glm::vec3(0.f, 0.5f, 1.f);
	incidentRay.direction = glm::vec3(0.f, 0.f, -1.f);
	incidentRay.direction = glm::normalize(incidentRay.direction);

	std::cout << " --- Incident ray onto glass sphere from outside ---" << std::endl;
	{
		HitRecord hitRec;
		bool hit = s.intersect(incidentRay, hitRec, "primary", true);
		assert(true == hit);
		assert(2 == hitRec.numHits);
		assert(true == hitRec.front_face);

		//Calculate fresnel number at hitpoint
		// We have our incident ray and its normalized
		// We have the normal at the hitpoint and its normalized
		// We know the ray is in medium with a low refracive index entering a medium with a higher refractive index
		// So we treat this as entering from a new medium
		float cosi = glm::clamp(-1.f, 1.f, glm::dot(incidentRay.direction, hitRec.hitNormal));
		assert(cosi < 0.f);

		float etai = 1.0f;
        float etat = s.IndexOfRefraction;

		float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
		bool totalInternalReflection = sint >= 1.f ? true : false;
		// If sint is greater than or equal to 1, then we have total internal reflection
		// Given this test incident ray and the angle at which it hits the sphere, it should NOT be the case
		assert(false == totalInternalReflection);

		//Calculate the frenel number
		float cost = sqrtf(std::max(0.f, 1.f - sint * sint)); 
		cosi = fabsf(cosi); 
		float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost)); 
		float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost)); 
		float kr = 0.5f * (Rs * Rs + Rp * Rp);

		std::cout << "Kr value is " << kr << std::endl << std::flush;

		// Check that the reflectance coefficient is between 0 and 1.
		// its possbile the ray fully refracts. In wich case kr = 0
		assert(kr > 0.f && kr < 1.f);

		// Check this kr value against actual fresnel function call
		Scene scene;
		float fresnelKr;
		scene.fresnel(incidentRay.direction, hitRec.hitNormal, s.IndexOfRefraction, fresnelKr, true);
		assert( kr == fresnelKr );
	}
	std::cout << " --- Refracted ray total internal reflection ---" << std::endl;
	{
		HitRecord hitRec;
		bool hit = s.intersect(incidentRay, hitRec, "primary", true);

		glm::vec3 bias = 0.00001f * hitRec.hitNormal;

		// Becuase our refracted ray is entering from outside into a sphere,
		// we push the origin of the refracted ray just inside the sphere with a small bias
		Ray refractedRay;
		refractedRay.origin = hitRec.hitPoint - bias; // bias is in direction of normal, so we subract to get inside the sphere
		refractedRay.direction = glm::refract(incidentRay.direction, hitRec.hitNormal, s.IndexOfRefraction);
		std::cout << "Incident ray  = " << printVec(incidentRay.direction) << std::endl << std::flush;
		std::cout << "Normal        = " << printVec(hitRec.hitNormal) << std::endl << std::flush;
		std::cout << "Refracted ray = " << printVec(refractedRay.direction) << std::endl << std::flush;
		std::cout << "Ref ray magnitude = " << glm::length(refractedRay.direction) << std::endl << std::flush;

		hitRec.front_face = true;
		hitRec.numHits = 0;
		hitRec.t = numeric_limits<float>::infinity();
		hit = s.intersect(refractedRay, hitRec, "Rr", true);
		assert(true == hit);
		assert(2 == hitRec.numHits);
		assert(false == hitRec.front_face);

		//Calculate fresnel number at hitpoint
		// We have our incident ray and its normalized
		// We have the normal at the hitpoint and its normalized
		// We know the ray is in medium with a high refracive index entering a medium with a lower refractive index
		// So we treat this as entering from an existing medium to air
		float cosi = glm::clamp(-1.f, 1.f, glm::dot(refractedRay.direction, hitRec.hitNormal));
		
		float etai = 1.0f;
        float etat = s.IndexOfRefraction;

		assert(cosi >= 0.f); // We are inside the high refractive index object
		std::swap(etai, etat); // so swap the etas
		
		float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
		bool totalInternalReflection = sint >= 1.f ? true : false;
		// If sint is greater than or equal to 1, then we have total internal reflection
		// Given this test incident ray and the angle at which it hits the sphere, it should NOT be the case
		assert(true == totalInternalReflection);

		float kr = 1.f;
		std::cout << "Kr value is " << kr << std::endl << std::flush;

		// Check this kr value against actual fresnel function call
		Scene scene;
		float fresnelKr;
		scene.fresnel(refractedRay.direction, hitRec.hitNormal, s.IndexOfRefraction, fresnelKr, true);
		assert( kr == fresnelKr );
	}
	std::cout << " --- END FRESNEL TEST ---" << std::endl;
}

void TestShallowAngleRefractedRay()
{
	float eps = std::numeric_limits<float>::epsilon();
	float eps2 = eps * 2.0f;
	std::cout << " --- BEGIN SHALLOW ANGLE RAY TEST ---" << std::endl;

	Scene scene;
	// Create a sphere with refractive index near that of glass
	Sphere s;
	s.IndexOfRefraction = 1.5;

	// Create a ray that will intersect with the sphere at a 45 degree angle
	Ray incidentRay;
	incidentRay.origin = glm::vec3(0.f, 0.95f, 1.f);
	incidentRay.direction = glm::vec3(0.f, 0.f, -1.f);
	incidentRay.direction = glm::normalize(incidentRay.direction);

	Ray previousRay;
	std::cout << " --- Shallow Angle Refracted ray -> reflected ray ---" << std::endl << std::flush;
	{
		HitRecord hitRec;
		bool hit = s.intersect(incidentRay, hitRec, "primary", false);
		assert(true == hit);
		assert(2 == hitRec.numHits);
		assert(true == hitRec.front_face);
		previousRay = incidentRay;
		glm::vec3 bias = 0.00001f * hitRec.hitNormal;

		float kr = 0.f;
		scene.fresnel(incidentRay.direction, hitRec.hitNormal, s.IndexOfRefraction, kr, false);
		std::cout << "Kr value is " << kr << std::endl << std::flush;

		// Becuase our refracted ray is entering from outside into a sphere,
		// we push the origin of the refracted ray just inside the sphere with a small bias
		Ray refractedRay;
		refractedRay.origin = hitRec.hitPoint - bias; // bias is in direction of normal, so we subract to get inside the sphere
		refractedRay.direction = glm::refract(incidentRay.direction, hitRec.hitNormal, s.IndexOfRefraction);
		std::cout << "refracted Ray  = " << printVec(refractedRay.direction) << std::endl << std::flush;

		// Test that refracted ray is valid
		if(glm::length(refractedRay.direction) >= 1.f - eps)
		{
			// reset hitRec
			hitRec.hitNormal = glm::vec3(0.f);
			hitRec.hitPoint = glm::vec3(0.f);
			hitRec.front_face = true;
			hitRec.numHits = 0;
			hitRec.t = numeric_limits<float>::infinity();
			hit = s.intersect(refractedRay, hitRec, "Rr", false);

			kr = 0.f;
			scene.fresnel(refractedRay.direction, hitRec.hitNormal, s.IndexOfRefraction, kr, false);
			std::cout << "Kr value is " << kr << std::endl << std::flush;
			previousRay = refractedRay;
		}
		else
		{
			std::cout << "refracted Ray is not a valid ray" << std::endl << std::flush;
		}

		int reflectedRayCount = 0;
		//Reflect the ray
		Ray reflectedRay;
		while(kr > 0.05f)
		{
			reflectedRayCount++;
			//std::cout << "Reflected Ray " << reflectedRayCount << std::endl << std::flush;
			bias = 0.00001f * hitRec.hitNormal;
			
			reflectedRay.origin = hitRec.front_face ? hitRec.hitPoint + bias : hitRec.hitPoint - bias;
			reflectedRay.direction = glm::reflect(previousRay.direction, hitRec.front_face ? hitRec.hitNormal: -hitRec.hitNormal);
			std::cout << "Normal        = " << printVec(hitRec.hitNormal) << std::endl << std::flush;
			std::cout << "RefLected ray = " << printVec(reflectedRay.direction) << std::endl << std::flush;
			std::cout << "RefLected ray magnitude = " << glm::length(reflectedRay.direction) << std::endl << std::flush;
			
			if(glm::length(reflectedRay.direction) > 1.f - eps)
			{
				// reset hitRec
				hitRec.hitNormal = glm::vec3(0.f);
				hitRec.hitPoint = glm::vec3(0.f);
				hitRec.front_face = true;
				hitRec.numHits = 0;
				hitRec.t = numeric_limits<float>::infinity();

				// Test intersection
				hit = s.intersect(reflectedRay, hitRec, "Rl", true);
				
				// calculate fresnel
				kr = 0.f;
				scene.fresnel(reflectedRay.direction, hitRec.hitNormal, s.IndexOfRefraction, kr, false);
				std::cout << "Kr value is " << kr << std::endl << std::flush;
				previousRay = reflectedRay;
			}
			else
			{
				std::cout << "reflected ray is not a valid ray " << std::endl << std::flush;
				break;
			}
		}
		std::cout << "Num reflected Rays = " << reflectedRayCount << std::endl << std::flush;
		// calculate new ray(s)
		if(kr < 1.0f) // we have to refract
		{
			//refract
		}
		//reflect
	}
	std::cout << " --- END SHALLOW ANGLE RAY TEST ---" << std::endl;
}

void Test45DegreeRefractedRay()
{
	float eps = std::numeric_limits<float>::epsilon();
	float eps2 = eps * 2.0f;
	std::cout << " --- BEGIN 45 DEGREE RAY TEST ---" << std::endl;

	Scene scene;
	// Create a sphere with refractive index near that of glass
	Sphere s;
	s.IndexOfRefraction = 1.5;

	// Create a ray that will intersect with the sphere at a 45 degree angle
	Ray incidentRay;
	incidentRay.origin = glm::vec3(0.f, 0.5f * sqrt(2.f) - 0.250f, 1.f);
	incidentRay.direction = glm::vec3(0.f, 0.f, -1.f);
	incidentRay.direction = glm::normalize(incidentRay.direction);

	Ray previousRay;
	std::cout << " --- 45 DEGREE  Refracted ray -> reflected ray ---" << std::endl << std::flush;
	{
		HitRecord hitRec;
		bool hit = s.intersect(incidentRay, hitRec, "primary", false);
		assert(true == hit);
		assert(2 == hitRec.numHits);
		assert(true == hitRec.front_face);
		previousRay = incidentRay;

		glm::vec3 bias = 0.00001f * hitRec.hitNormal;

		float kr = 0.f;
		scene.fresnel(incidentRay.direction, hitRec.hitNormal, s.IndexOfRefraction, kr, false);
		std::cout << "Kr value is " << kr << std::endl << std::flush;

		// Becuase our refracted ray is entering from outside into a sphere,
		// we push the origin of the refracted ray just inside the sphere with a small bias
		Ray refractedRay;
		refractedRay.origin = hitRec.hitPoint - bias; // bias is in direction of normal, so we subract to get inside the sphere
		refractedRay.direction = glm::refract(previousRay.direction, hitRec.hitNormal, s.IndexOfRefraction);
		std::cout << "refracted Ray  = " << printVec(refractedRay.direction) << std::endl << std::flush;

		// Test that refracted ray is valid
		if(glm::length(refractedRay.direction) >= 1.f - eps)
		{
			// reset hitRec
			hitRec.hitNormal = glm::vec3(0.f);
			hitRec.hitPoint = glm::vec3(0.f);
			hitRec.front_face = true;
			hitRec.numHits = 0;
			hitRec.t = numeric_limits<float>::infinity();
			hit = s.intersect(refractedRay, hitRec, "Rr", false);

			kr = 0.f;
			scene.fresnel(refractedRay.direction, hitRec.hitNormal, s.IndexOfRefraction, kr, false);
			std::cout << "Kr value is " << kr << std::endl << std::flush;
			previousRay = refractedRay;
		}
		else
		{
			std::cout << "refracted Ray is not a valid ray" << std::endl << std::flush;
		}

		int reflectedRayCount = 0;
		//Reflect the ray
		Ray reflectedRay;
		while(kr > 0.05f)
		{
			reflectedRayCount++;
			//std::cout << "Reflected Ray " << reflectedRayCount << std::endl << std::flush;
			bias = 0.00001f * hitRec.hitNormal;
			
			reflectedRay.origin = hitRec.front_face ? hitRec.hitPoint + bias : hitRec.hitPoint - bias;
			reflectedRay.direction = glm::reflect(previousRay.direction, hitRec.front_face ? hitRec.hitNormal: -hitRec.hitNormal);
			std::cout << "Normal        = " << printVec(hitRec.hitNormal) << std::endl << std::flush;
			std::cout << "RefLected ray = " << printVec(reflectedRay.direction) << std::endl << std::flush;
			std::cout << "RefLected ray magnitude = " << glm::length(reflectedRay.direction) << std::endl << std::flush;
			
			if(glm::length(reflectedRay.direction) > 1.f - eps)
			{
				// reset hitRec
				hitRec.hitNormal = glm::vec3(0.f);
				hitRec.hitPoint = glm::vec3(0.f);
				hitRec.front_face = true;
				hitRec.numHits = 0;
				hitRec.t = numeric_limits<float>::infinity();

				// Test intersection
				hit = s.intersect(reflectedRay, hitRec, "Rl", true);
				
				// calculate fresnel
				kr = 0.f;
				scene.fresnel(reflectedRay.direction, hitRec.hitNormal, s.IndexOfRefraction, kr, false);
				std::cout << "Kr value is " << kr << std::endl << std::flush;
				previousRay = reflectedRay;
			}
			else
			{
				std::cout << "reflected ray is not a valid ray " << std::endl << std::flush;
				break;
			}
		}
		std::cout << "Num reflected Rays = " << reflectedRayCount << std::endl << std::flush;
		// calculate new ray(s)
		if(kr < 1.0f) // we have to refract
		{
			refractedRay.origin = hitRec.hitPoint - bias; // bias is in direction of normal, so we subract to get inside the sphere
			refractedRay.direction = glm::refract(incidentRay.direction, hitRec.hitNormal, s.IndexOfRefraction);
			std::cout << "refracted Ray  = " << printVec(refractedRay.direction) << std::endl << std::flush;
		}
		//reflect
	}
	std::cout << " --- END 45 DEGREE RAY TEST ---" << std::endl;
}

void TestSteepAngleRefractedRay()
{
	float eps = std::numeric_limits<float>::epsilon();
	float eps2 = eps * 2.0f;
	std::cout << " --- BEGIN STEEP ANGLE RAY TEST ---" << std::endl;

	Scene scene;
	// Create a sphere with refractive index near that of glass
	Sphere s;
	s.IndexOfRefraction = 1.5;

	// Create a ray that will intersect with the sphere at a 45 degree angle
	Ray incidentRay;
	incidentRay.origin = glm::vec3(0.f, 0.05f, 1.f);
	incidentRay.direction = glm::vec3(0.f, 0.f, -1.f);
	incidentRay.direction = glm::normalize(incidentRay.direction);

	Ray previousRay;
	int reflectedRayCount = 0;
	int refractedRayCount = 0;
	std::cout << " --- Steep Angle  Refracted ray -> reflected ray ---" << std::endl << std::flush;
	{
		HitRecord hitRec;
		bool hit = s.intersect(incidentRay, hitRec, "primary", true);
		assert(true == hit);
		assert(2 == hitRec.numHits);
		assert(true == hitRec.front_face);
		previousRay = incidentRay;
		glm::vec3 bias = 0.00001f * hitRec.hitNormal;

		float kr = 0.f;
		scene.fresnel(incidentRay.direction, hitRec.hitNormal, s.IndexOfRefraction, kr, false);
		std::cout << "Kr value is " << kr << std::endl << std::flush;

		// Becuase our refracted ray is entering from outside into a sphere,
		// we push the origin of the refracted ray just inside the sphere with a small bias
		Ray refractedRay;
		refractedRay.origin = hitRec.front_face ? hitRec.hitPoint - bias : hitRec.hitPoint + bias; // bias is in direction of normal, so we subract to get inside the sphere
		refractedRay.direction = glm::refract(previousRay.direction, hitRec.hitNormal, s.IndexOfRefraction);
		std::cout << "refracted Ray  = " << printVec(refractedRay.direction) << std::endl << std::flush;

		// Test that refracted ray is valid
		if(glm::length(refractedRay.direction) >= 1.f - eps)
		{
			refractedRayCount++;
			// reset hitRec
			hitRec.hitNormal = glm::vec3(0.f);
			hitRec.hitPoint = glm::vec3(0.f);
			hitRec.front_face = true;
			hitRec.numHits = 0;
			hitRec.t = numeric_limits<float>::infinity();
			hit = s.intersect(refractedRay, hitRec, "Rr", true);

			kr = 0.f;
			scene.fresnel(refractedRay.direction, hitRec.hitNormal, s.IndexOfRefraction, kr, false);
			std::cout << "Kr value is " << kr << std::endl << std::flush;
			previousRay = refractedRay;
		}
		else
		{
			std::cout << "refracted Ray is not a valid ray" << std::endl << std::flush;
		}

		//Reflect the ray
		Ray reflectedRay;
		while(kr > 0.05f)
		{
			//std::cout << "Reflected Ray " << reflectedRayCount << std::endl << std::flush;
			bias = 0.00001f * hitRec.hitNormal;
			
			reflectedRay.origin = hitRec.front_face ? hitRec.hitPoint + bias : hitRec.hitPoint - bias;
			reflectedRay.direction = glm::reflect(previousRay.direction, hitRec.front_face ? hitRec.hitNormal: -hitRec.hitNormal);
			std::cout << "Normal        = " << printVec(hitRec.hitNormal) << std::endl << std::flush;
			std::cout << "RefLected ray = " << printVec(reflectedRay.direction) << std::endl << std::flush;
			std::cout << "RefLected ray magnitude = " << glm::length(reflectedRay.direction) << std::endl << std::flush;
			
			if(glm::length(reflectedRay.direction) > 1.f - eps)
			{
				reflectedRayCount++;
				// reset hitRec
				hitRec.hitNormal = glm::vec3(0.f);
				hitRec.hitPoint = glm::vec3(0.f);
				hitRec.front_face = true;
				hitRec.numHits = 0;
				hitRec.t = numeric_limits<float>::infinity();

				// Test intersection
				hit = s.intersect(reflectedRay, hitRec, "Rl", true);
				
				// calculate fresnel
				kr = 0.f;
				scene.fresnel(reflectedRay.direction, hitRec.hitNormal, s.IndexOfRefraction, kr, false);
				std::cout << "Kr value is " << kr << std::endl << std::flush;
				previousRay = reflectedRay;
			}
			else
			{
				std::cout << "reflected ray is not a valid ray " << std::endl << std::flush;
				break;
			}
		}
		std::cout << "Num reflected Rays = " << reflectedRayCount << std::endl << std::flush;
		// calculate new ray(s)
		if(kr < 1.0f) // we have to refract
		{
			glm::vec3 bias = 0.00001f * hitRec.hitNormal;

			refractedRay.origin = hitRec.front_face ? hitRec.hitPoint - bias : hitRec.hitPoint + bias; // bias is in direction of normal, so we subract to get inside the sphere
			refractedRay.direction = glm::refract(previousRay.direction, -hitRec.hitNormal, s.IndexOfRefraction);
			std::cout << "refracted Ray  = " << printVec(refractedRay.direction) << std::endl << std::flush;

			// Test that refracted ray is valid
			if(glm::length(refractedRay.direction) >= 1.f - eps)
			{
				refractedRayCount++;
				// reset hitRec
				hitRec.hitNormal = glm::vec3(0.f);
				hitRec.hitPoint = glm::vec3(0.f);
				hitRec.front_face = true;
				hitRec.numHits = 0;
				hitRec.t = numeric_limits<float>::infinity();
				hit = s.intersect(refractedRay, hitRec, "Rr", true);

				if(hitRec.numHits > 0)
				{
					kr = 0.f;
					scene.fresnel(refractedRay.direction, hitRec.hitNormal, s.IndexOfRefraction, kr, false);
					std::cout << "Kr value is " << kr << std::endl << std::flush;
					previousRay = refractedRay;
				}
			}
			else
			{
				std::cout << "refracted Ray is not a valid ray" << std::endl << std::flush;
			}
		}
		//reflect
	}
	std::cout << " --- END STEEP ANGLE RAY TEST ---" << std::endl;
}

void TestIndexOfRefraction()
{

}

int main()
{
	//TestRaySphereIntersection();
	//TestReflectDirection();

	//TestFresnelNumber();
	//TestShallowAngleRefractedRay();
	//Test45DegreeRefractedRay();
	//TestSteepAngleRefractedRay();
	// TestIndexOfRefraction();

	Scene scene;
	scene.init(1.5f);

	Renderer r;
	r.init();
	r.setFramebufferDimensions(1920, 1080);
	scene.rayTracer = r;

	scene.render();
	/*for( float ior = 1.1f; ior < 36.0f; ior += 0.5f)
	{
		scene.clearScene();
		scene.init(ior);
		scene.render();
	}*/
	return 0;
}