#pragma once

#include <iostream>
#include "ObjectBase.h"
#include <glm/glm.hpp>

using namespace std;

struct HitRecord
{
    /*
    glm::vec3 RayOrigin;
    glm::vec3 rayDirection;
    uint32_t intersectionCount;
    float * intersections_t;
    ObjectBase* intersectedObject;
    */

   glm::vec3 hitPoint;
   glm::vec3 hitNormal;
   uint32_t numHits;
   float t;
   bool front_face;

   std::string printVec3(const glm::vec3& vec) const
   {
       std::stringstream ss;
       ss << "(" << vec.x << "," << vec.y << "," << vec.z << ")";
       return ss.str();
   }

   void print() const
   {
       std::cout << "hitPoint = " << printVec3(hitPoint) << std::endl;
       std::cout << "hitNormal = " << printVec3(hitNormal) << std::endl;
       std::cout << "Normal Length = " << glm::length(hitNormal) << std::endl;
       std::cout << "time t = " << t << std::endl;
       std::cout << "front_face = " << front_face << std::endl;
   }
};