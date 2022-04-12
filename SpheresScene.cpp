#include "SpheresScene.h"

#include <cassert>

using namespace std;

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

	{
		HitRecord rec;
		bool hit = s.intersect(ray, rec, false);
		assert(hit != true);
	}

	// Now create a ray this is exactly +1  on the y axis
	// This ray should intersect the sphere tangentailly
	ray.origin = glm::vec3(0.f, 1.f, -5.0f);
	{
		HitRecord rec;
		bool hit = s.intersect(ray, rec, false);
		assert(hit == true);
		assert(rec.numHits == 1);
		assert(rec.front_face == true);

		// Check that the normal has been normalized
		assert(glm::length(rec.hitNormal) >= 0.9999f);
		//Check that the ray direction and the sphere normal are negative
		assert(glm::dot(ray.direction, rec.hitNormal) <= 0.f);
	}

	// Now create a ray that will interserct the sphere twice
	ray.origin = glm::vec3(0.f, 0.0f, -5.0f);
	{
		HitRecord rec;
		bool hit = s.intersect(ray, rec, true);
		assert(hit == true);
		assert(rec.numHits == 2);
		assert(rec.front_face == true);

		// Check that the normal has been normalized
		assert(glm::length(rec.hitNormal) >= 0.9999f);
		//Check that the ray direction and the sphere normal are negative
		assert(glm::dot(ray.direction, rec.hitNormal) <= 0.f);
	}

	// Now create a ray whos' origin is inside the sphere
	ray.origin = glm::vec3(0.f, 0.f, 0.f);
	{
		HitRecord rec;
		bool hit = s.intersect(ray, rec, true);
		assert(hit == true);
		assert(rec.numHits == 2);
		assert(rec.front_face == false);

		// Check that the normal has been normalized
		assert(glm::length(rec.hitNormal) >= 0.9999f);
		//Check that the ray direction and the sphere normal are negative
		assert(glm::dot(ray.direction, rec.hitNormal) <= 0.f);
	}
}

int main()
{
	TestRaySphereIntersection();

	Scene scene;
	scene.init(1.5f);

	Renderer r;
	r.init();
	r.setFramebufferDimensions(640, 480);
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