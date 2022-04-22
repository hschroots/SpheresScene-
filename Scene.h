#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include <exception>
#include <vector>
#include <numeric>
#include <random>

#include "Renderer.h"
#include "Camera.h"
#include "ObjectBase.h"
#include "Sphere.h"
#include "Plane.h"
#include "Light.h"
#include "Ray.h"

#define MAX_TRACE_DEPTH 32
#define RAYS_PER_PIXEL 8

uint32_t primaryCount = 0;
uint32_t reflectCount = 0;
uint32_t refractCount = 0;

glm::vec3 clearColor(0.3f, 0.9f, 1.0f);

class Scene
{
    public:
    Scene() : initialized(false) { }
    ~Scene() {}

    void init(float ior = 1.1f)
    {
        if(!initialized)
        {
            setupLights();
            setupGeometry(ior);
            initialized = true;
        }
    }
    void clearScene()
    {
        lights.clear();
        for( int i = 0; i < objects.size(); i++)
        {
            delete objects[i];
        }
        objects.clear();
        initialized = false;
    }

    inline double random_double() {
        static std::uniform_real_distribution<double> distribution(0.0, 1.0);
        static std::mt19937 generator;
        return distribution(generator);
    }

    void render()
    {
        if(initialized)
        {

        if(rayTracer.width < 1)
        {
            std::cerr << "Renderer not set" << std::endl;
            return;
        }

        uint32_t width = rayTracer.width;
        uint32_t height = rayTracer.height;

        // perterb ray within pixel boundaries
        // run 4 rays per pixel
        // ray offset to lower left
        float dx = 1.f / (float)(width-1);
        float dy = 1.f / (float)(height - 1);


        offsets[0] = glm::vec3(-0.166667f * dx, -0.166667f * dy, 0.0f);
        offsets[1] = glm::vec3(           0.0f,  0.333333f * dy, 0.0f);
        offsets[2] = glm::vec3( 0.416667f * dx, 0.0833333f * dy, 0.0f);
        offsets[3] = glm::vec3(           0.0f,       0.5f * dy, 0.0f);

        Ray ray;
        for(int row = 0; row < height; row++)
            for(int col = 0; col < width; col++)
            {
                bool debug = false;
                glm::vec3 colorVal = glm::vec3(0);

                //for(int s = 0; s < RAYS_PER_PIXEL; s++)
                float u;
                float v;
                int sqrtfloor = std::max(1.f,(float)std::floor(std::sqrt(RAYS_PER_PIXEL)));
                float delta = 1.f / (float)sqrtfloor;
                for(int s = 0; s < sqrtfloor; s++)
                for(int t = 0; t < sqrtfloor; t++)
                {
                    if(RAYS_PER_PIXEL > 1)
                    {
                        //u = (col + random_double()) / (width-1);
                        //v = (row + random_double()) / (height-1);
                        u = (col + delta * t) / (width-1);
                        v = (row + delta * s) / (height-1);
                    }else{
                        u = (float)(col) / (float)(width-1);
                        v = (float)(row) / (float)(height-1);
                    }
                    rayTracer.computeRay(u, v, ray);

                    // refractions accurs @ row 240 from col 334 to 549
                    // particularly col 
                    if( col == 500 && row == 240)
                    {
                        debug = true;
                    }
                    if(debug)
                    {
                        std::cout << "(row,col, r) = (" << row << "," << col << "," << s << ")" << std::endl;
                        std::cout << "RayOrigin: " << ray.origin.x << ", " << ray.origin.y << ", " << ray.origin.z << std::endl;
                        std::cout << "Ray: " << ray.direction.x << ", " << ray.direction.y << ", " << ray.direction.z << std::endl;
                    }
                    primaryCount++;
                    colorVal += traceRay(ray, 0, col, row, "primary", 1.0f, debug);
                }
                float scale = 1.f/(float)RAYS_PER_PIXEL;
                colorVal = scale * colorVal;

                /*
                float power = 1.f/2.4f;
                colorVal.b = colorVal.b <= 0.0031308f ? 12.92 * colorVal.b : 1.055f * std::pow(colorVal.b, power) - 0.055;
                colorVal.g = colorVal.g <= 0.0031308f ? 12.92 * colorVal.g : 1.055f * std::pow(colorVal.g, power) - 0.055;
                colorVal.r = colorVal.r <= 0.0031308f ? 12.92 * colorVal.r : 1.055f * std::pow(colorVal.r, power) - 0.055;
                */
                // Gamma correct
                /*
                colorVal.b = sqrt(colorVal.b);
                colorVal.g = sqrt(colorVal.g);
                colorVal.r = sqrt(colorVal.r);
                */

                colorVal = glm::clamp(colorVal, 0.f, 1.f);

                const float b = colorVal.b;
                const float g = colorVal.g;
                const float r = colorVal.r;

                //if(col % (width -1) == 0)
                    //std::cout << "Writing pixel (" << col << "," << row << ")" << std::endl << std::flush;
                rayTracer.framebuffer.SetPixel(Pixel(b,g,r), col, row);
            }

        std::cout << "Counters: " << std::endl;
        std::cout << "Primary Rays  = " << primaryCount << std::endl;
        std::cout << "Refract Rays  = " << refractCount << std::endl;
        std::cout << "Reflect Rays  = " << reflectCount << std::endl;
        std::stringstream filename;
        filename << "out/images/output_ior_" << objects[0]->IndexOfRefraction << "_" << rayTracer.width << "x" << rayTracer.height << ".bmp";
        rayTracer.framebuffer.Export(filename.str().c_str());
        }
    }

    float FresnelReflectAmount(float n1, float n2, glm::vec3 normal, glm::vec3 incident, float reflectivity)
    { 
        // https://blog.demofox.org/2017/01/09/raytracing-reflection-refraction-fresnel-total-internal-reflection-and-beers-law/
        // Schlick aproximation
        float r0 = (n1-n2) / (n1+n2);
        r0 *= r0;
        float cosX = -glm::dot(normal, incident);
        if (n1 > n2)
        {
            float n = n1/n2;
            float sinT2 = n*n*(1.0-cosX*cosX);
            // Total internal reflection
            if (sinT2 > 1.0)
                return 1.0;
            cosX = sqrt(1.0-sinT2);
        }
        float x = 1.0-cosX;
        float ret = r0+(1.0-r0)*x*x*x*x*x;

        // adjust reflect multiplier for object reflectivity
        ret = (reflectivity + (1.0-reflectivity) * ret);
        return ret;
    }

    void mySchlick(glm::vec3 const I, glm::vec3 const N, bool inside, float ior, float& kr)
    {
        float cosi = glm::dot(N,I);

        float eta = 1.f / ior;
        if(inside)
        {
            eta = ior;
            float sin2_t = eta * eta * ( 1.f - cosi * cosi);
            if(sin2_t > 1.0f)
            {
                kr = 1.f;
                return;
            }
            float cost = std::sqrt(1.f - sin2_t);
            cosi = cost;
        }

        float c = (1.f - ior) / (1.f + ior);
        float r0 = c * c;

        float q = (1.f - cosi);
        float q5 = q * q * q * q * q;
        kr = r0 + ( 1 - r0 ) * q5;
    }
    void fresnel(glm::vec3& incidentRay, glm::vec3& norm, float& ior, float& kr, bool debug = false)
    {
        glm::vec3 I = glm::normalize(incidentRay);
        glm::vec3 N = glm::normalize(norm);

        float cosi = glm::clamp(-1.0f, 1.0f, glm::dot(N, I));
        float etai = 1.0f;
        float etat = ior;
        
        if(cosi > 0.0f) // we are inside the object
        {
            std::swap(etai, etat);
        }
        float sint = etai / etat * std::sqrt(std::max(0.f, 1.f - cosi * cosi)); 
        // Total internal reflection
        if (sint >= 1.f) { 
            kr = 1.f;
        } 
        else { 
            float cost = std::sqrt(std::max(0.f, 1.f - sint * sint)); 
            cosi = fabsf(cosi); 
            float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost)); 
            float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost)); 
            kr = 0.5f * (Rs * Rs + Rp * Rp); 
        }
    }

    glm::vec3 refract(glm::vec3 incidentRay, glm::vec3 norm, float ior, bool debug = false)
    {
        glm::vec3 I = glm::normalize(incidentRay);
        glm::vec3 N = glm::normalize(norm);

        float cosi = glm::clamp(-1.0f, 1.0f, glm::dot(N, I));
        float etai = 1.0f;
        float etat = ior;

        if(cosi < 0.f)
        {
            cosi = -cosi;
        }else{
            //reverse the normal because we are inside the object
            N = -N;
            //swap the index of refractions
            std::swap(etai, etat);
        }
        float eta = etai / etat;
        //Check to see if we have total internal reflection
        float k = 1.0f - eta * eta * (1.0f - cosi * cosi);
        if( k < 0.f )
        {
            // total internal refelection
            // Return white light
            return glm::vec3(1.0f);
        }
        else
        {
            return eta * I - ( eta * cosi + std::sqrt(k) ) * N;
        }
    }

    glm::vec3 reflect(glm::vec3 incidentRay, glm::vec3 norm, bool debug = false)
    {
        glm::vec3 I = glm::normalize(incidentRay);
        glm::vec3 N = glm::normalize(norm);
        return I - N * glm::dot(N, I) * 2.0f;
    }

    float traceShadowRay(Ray& ray, bool debug)
    {
        for(int idx = 0; idx < objects.size(); idx++)
        {
            HitRecord rec;
            if(rayTracer.intersect(ray, objects[idx], rec, "shadow", debug))
            {
                return 1.0f;
            }
        }
        return 0.0f;
    }

    glm::vec3 traceRay(Ray& ray, uint32_t traceDepth, int col, int row, string type, float current_ior, bool debug)
    {
        if(traceDepth > MAX_TRACE_DEPTH)
            return clearColor; //glm::vec3(1.0f);

        glm::vec3 color(0.0f);
        ObjectBase* obj = 0;
        float minDist = std::numeric_limits<float>::max();

        HitRecord hitRec;
        for(int idx = 0; idx < objects.size(); idx++)
        {
            HitRecord rec;
            if(rayTracer.intersect(ray, objects[idx], rec, type, debug))
            {
                float distance = glm::distance(ray.origin, rec.hitPoint);
                if(debug)
                {
                    std::cout << type + " " << "Intersected with object idx =  " << idx << std::endl;
                    std::cout << type + " "<< "Object is at a distance of  " << distance << std::endl;
                    std::cout << type + " "<< "Is that less than the current distance of " << minDist << std::endl;
                }
                if(distance < minDist)
                {
                    obj = objects[idx];
                    hitRec = rec;
                    if (debug && glm::length(hitRec.hitNormal) < 0.99999f)
                    {
                        std::cout << type + " "<< "Pixel (col,row) =  " << col << "," << row << std::endl;
                        
                        hitRec.print();
                        std::cout << type + " "<< "Object Idx = " << idx << std::endl;
                    }
                    minDist = distance;
                }
            }
        }

        // If the min distance is still some large number, then the ray didn't hit anything
        // Set the color to the background color
        if(minDist >= 9999999999.99f)
        {
            if(debug)
            {
                std::cout << type + " " << "Missed all objects. Returning sky color " << std::endl;
            }
            return clearColor;
        }
        else // we hit an object and it is the closest object
        {
            if(col == 388 && row == 240 && type == "primary")
            {
                std::cout << type + " " << "Debug value is " << debug << std::endl << std::flush;
                debug = true;
            }
            glm::vec3 bias = 0.00001f * hitRec.hitNormal;

            /// ---------------- SHADOWS --------------------- ///
            Light light = lights[0];
            glm::vec3 lightdir = light.position  - hitRec.hitPoint;
            lightdir = glm::normalize(lightdir);

            float shadowK = 0.6f;
            float shadowAcc = 0.f;
            int numShadowRays = 4;
            float invShadowRays = 1.f/(float)numShadowRays;
            
            // Only shadow the floor
            if (dynamic_cast<Plane*>(obj) != nullptr)
            {
                if(type == "primary" && hitRec.front_face)
                {
                    for( int i = 0; i < numShadowRays; i++)
                    {
                        Ray shadowRay;
                        shadowRay.origin = hitRec.hitPoint + bias;
                        glm::vec3 shadowToLight = light.position  - shadowRay.origin;
                        float distanceToLight = glm::length(shadowToLight);
                        shadowToLight = glm::normalize(shadowToLight);
                        shadowRay.direction = shadowToLight; //glm::normalize(lightdir + offsets[i % 4]);
                        shadowRay.tmin = 2.0f * 0.00001f;
                        shadowRay.tmax = distanceToLight - 0.00002f;
                        shadowAcc += invShadowRays * traceShadowRay(shadowRay, debug);
                    }
                }
            }
            shadowK *= shadowAcc;

            /// ------------  FRESNEL REFLECTIVITY --------------- ///
            //Calculate reflectivity
            float kr = 0.f;
            float current_ior = hitRec.front_face ? 1.f : obj->IndexOfRefraction;
            float next_ior = hitRec.front_face ? obj->IndexOfRefraction : 1.f;
            float reflectiveK = FresnelReflectAmount(current_ior, next_ior, hitRec.hitNormal, ray.direction, obj->reflectivity);

            /// ------------  CALCULATE ALBEDO -------------------- ////
            switch(obj ? obj->shadingType : -1)
            {
                case Opaque: // If the object is opaque shade it according to the blinn-phong model
                {
                    color = obj->ColorAt(hitRec.hitPoint);
                    // Apply blinn phong lighting model
                    color = BlinnPhong(color, light, hitRec, 100.f);
                }
                break;
                case Transparent:
                {
                    if(reflectiveK < 1.f - (2.f *numeric_limits<float>::epsilon()))
                    {
                        glm::vec3 refractColor = glm::vec3(1.f);
                        // Refract to get base color
                        Ray refractRay;
                        refractRay.origin = hitRec.front_face ? hitRec.hitPoint - bias : hitRec.hitPoint + bias; 
                        glm::vec3 local_norm = hitRec.front_face ? hitRec.hitNormal : -hitRec.hitNormal;
                        float eta = current_ior / next_ior; 
                        refractRay.direction = glm::refract(ray.direction, local_norm, eta);
                        refractRay.direction = glm::normalize(refractRay.direction);

                        bool refractRayValid = false;
                        if(glm::length(refractRay.direction) >= (1.f - (2.f * numeric_limits<float>::epsilon())))
                        {
                            refractCount++;
                            refractRayValid = true;
                            refractColor = traceRay(refractRay, traceDepth +1, col, row, "refract", 1.f, debug);
                        }
                        else{
                            std::cout << type + " " << "---- FAILED TO REFRACT RAY -------" << std::endl << std::flush;
                            glm::vec3 const N = local_norm;
                            glm::vec3 const I = ray.direction;
                            float const dotValue(dot(N, I));
                            float const k(static_cast<float>(1) - eta * eta * (static_cast<float>(1) - dotValue * dotValue));
                            bool const valid = k >= static_cast<float>(0);
                            std::cout << "Kr has value " << kr << " which means there exists a refracted ray, but refracted ray is " << printVec(refractRay.direction) << std::endl << std::flush;
                            std::cout << "Ray has length " << glm::length(refractRay.direction) << std::endl << std::flush;
                            std::cout << type + " " << "Current Ray direction = " << printVec(ray.direction) << std::endl << std::flush;
                            std::cout << type + " " << "hit normal = " << printVec(hitRec.hitNormal) << std::endl << std::flush;
                            std::cout << type + " " << "Front face is " << (hitRec.front_face ? "true" : "false") << std::endl << std::flush;
                            std::cout << type + " " << "Current ior = " << current_ior << std::endl << std::flush;
                            std::cout << type + " " << "Next ior = " << next_ior << std::endl << std::flush;
                            std::cout << type + " " << "Eta is " << eta << std::endl << std::flush;
                            std::cout << type + " " << "Dot value " << dotValue << std::endl << std::flush;
                            std::cout << type + " " << "k Value = " << k << std::endl << std::flush;
                            std::cout << type + " " << "k is " << (valid ? "valid" : "INVALID") << std::endl << std::flush;

                            //glm::vec3 const ReturnValue = valid ? (eta * I - (eta * dotValue + std::sqrt(k)) * N) : glm::vec3(0.f);
                        }
                        color = refractRayValid ? refractColor : obj->ColorAt(hitRec.hitPoint);
                        // Apply specular highlight
                        float specularK = SpecularK(lightdir, hitRec.hitPoint, hitRec.hitNormal, 100.f);
                        color = color + (specularK * light.color);  
                    }
                }
                break;
                /*
                case MyOldStuff: // If the object is transparent trace reflection and refraction rays
                {
                    glm::vec3 refractColor(1.f);
                    glm::vec3 reflectColor(1.f);
                    float kr = 0.0f;
                    float skr = 0.f;
                    // if fresnel number is 1, then we have total internal reflection
                    fresnel(ray.direction, hitRec.hitNormal, obj->IndexOfRefraction, kr, debug);
                    //schlick(ray.direction, hitRec.hitNormal, hitRec.front_face ? false : true, obj->IndexOfRefraction, skr);
                    //if( !(abs(kr - skr) < numeric_limits<float>::epsilon()) )
                    //{
                      //  std::cout << type + " " << "SCHLICK and FRESNEL NOT THE SAME" << std::endl;
                       // std::cout << type + " " << "SKR = " << skr << "  kr = " << kr << std::endl; 
                    //}
                    if(debug)
                    {
                        std::string status = (kr < 1.0f) ? " - refract " : " - total internal reflection ";
                        std::cout << type + " " << " kr = " << kr <<  status << std::endl << std::flush;
                        std::cout << type + " " << "ray direction = " << printVec(ray.direction) << std::endl << std::flush;
                        std::cout << type + " " << "hit normal = " << printVec(hitRec.hitNormal) << std::endl << std::flush;
                        std::cout << type + " " << "ray dot normal = " << glm::dot(ray.direction, hitRec.hitNormal) << std::endl << std::flush;
                    }
                    
                    if(kr < 1.f) // trace refraction
                    {
                        Ray refractRay;
                        // bias is in direction of normal, so we subract to get inside the sphere and add to get out of it
                        refractRay.origin = hitRec.front_face ? hitRec.hitPoint - bias : hitRec.hitPoint + bias; 
			            glm::vec3 local_norm = hitRec.front_face ? hitRec.hitNormal : -hitRec.hitNormal;
                        float eta = hitRec.front_face ? 1.f / obj->IndexOfRefraction : obj->IndexOfRefraction;
                        refractRay.direction = glm::refract(ray.direction, local_norm, eta);
                        refractRay.direction = glm::normalize(refractRay.direction);

                        // If the ray has length of almost 1
                        if(glm::length(refractRay.direction) >= (1.f - numeric_limits<float>::epsilon()))
                        {
                            if(debug)
                            {
                                //std::string status = (kr < 1.0f) ? " - refract " : " - total internal reflection ";
                                std::cout << type + " " << "Front face is " << (hitRec.front_face ? "true" : "false") << std::endl << std::flush;
                                std::cout << type + " " << "Refract Ray direction = " << printVec(refractRay.direction) << std::endl << std::flush;
                                //std::cout << type + " " << "hit normal = " << printVec(hitRec.hitNormal) << std::endl << std::flush;
                                //std::cout << type + " " << "ray dot normal = " << glm::dot(ray.direction, hitRec.hitNormal) << std::endl << std::flush;
                            }
                            //assert(glm::length(refractRay.direction) > 0.f);
                            //trace refracted ray
                            refractCount++;
                            refractColor = traceRay(refractRay, traceDepth +1, col, row, "refract", eta, debug);
                        }
                        else{
                            glm::vec3 const N = local_norm;
                            glm::vec3 const I = ray.direction;
                            float const dotValue(dot(N, I));
                            float const k(static_cast<float>(1) - eta * eta * (static_cast<float>(1) - dotValue * dotValue));
                            bool const valid = k >= static_cast<float>(0);
                            std::cout << "Kr has value " << kr << " which means there exists a refracted ray, but refracted ray is " << printVec(refractRay.direction) << std::endl << std::flush;
                            std::cout << "Ray has length " << glm::length(refractRay.direction) << std::endl << std::flush;
                            std::cout << type + " " << "Current Ray direction = " << printVec(ray.direction) << std::endl << std::flush;
                            std::cout << type + " " << "hit normal = " << printVec(hitRec.hitNormal) << std::endl << std::flush;
                            std::cout << type + " " << "Front face is " << (hitRec.front_face ? "true" : "false") << std::endl << std::flush;
                            std::cout << type + " " << " Eta is " << eta << std::endl << std::flush;
                            std::cout << type + " " << "Dot value " << dotValue << std::endl << std::flush;
                            std::cout << type + " " << "k Value = " << k << std::endl << std::flush;
                            std::cout << type + " " << "k is " << (valid ? "valid" : "INVALID") << std::endl << std::flush;

                            glm::vec3 const ReturnValue = valid ? (eta * I - (eta * dotValue + std::sqrt(k)) * N) : glm::vec3(0.f);
                        }   
                    }
                    if(kr > 0.05f)
                    {
                        //trace refelction
                        Ray reflectRay;
                        reflectRay.origin = hitRec.front_face ? hitRec.hitPoint + bias : hitRec.hitPoint - bias;
                        glm::vec3 nRef = hitRec.front_face ? hitRec.hitNormal : -hitRec.hitNormal;
                        reflectRay.direction = glm::reflect(ray.direction, nRef);

                        // If the ray has length of almost 1
                        if(glm::length(reflectRay.direction) >= (1.f - numeric_limits<float>::epsilon()))
                        {
                            reflectCount++;
                            reflectColor = traceRay(reflectRay, traceDepth + 1, col, row, "reflect", obj->IndexOfRefraction, debug);
                        }
                        else{
                            if(debug) std::cout << "Kr has value " << kr << " which means there exists a reflected ray, but reflected ray is " << printVec(reflectRay.direction) << std::endl << std::flush;
                        }
                    }
                    color += reflectColor * kr + refractColor * (1.f - kr);
                }
                break;
                */
                default:
                break;
            }
            if(hitRec.front_face)
            {
                if(reflectiveK > 0.f)
                {
                    glm::vec3 reflectColor(1.f);
                    Ray reflectRay;
                    reflectRay.origin = hitRec.front_face ? hitRec.hitPoint + bias : hitRec.hitPoint - bias;
                    glm::vec3 nRef = hitRec.front_face ? hitRec.hitNormal : -hitRec.hitNormal;
                    reflectRay.direction = glm::reflect(ray.direction, nRef);
                    reflectRay.direction = glm::normalize(reflectRay.direction);
                    reflectCount++;
                    reflectColor = traceRay(reflectRay, traceDepth + 1, col, row, "reflect", obj->IndexOfRefraction, debug);
                    reflectColor = BlinnPhong(reflectColor, light, hitRec, 100.f);
                    color = reflectiveK * reflectColor + (1.f-reflectiveK) * color;
                }
            }

            return (1.f-shadowK) * color;
        }
    }

    private:
    void setupLights()
    {
        Light light;
        light.position = glm::vec3(-1.f, 6.0f, -5.f);
        light.color = glm::vec3(1.0f, 1.0f, 1.0f);
        lights.push_back(light);
    }

    void setupGeometry(float ior = 1.1f)
    {
        glm::vec3 center;
        glm::vec3 scale;
        // Add sphere
        
        Sphere* s = new Sphere();
        center = glm::vec3(-2.0, 0.0f, -5.0f);
        s->transform.Translate(center);
        //glm::vec3 scale = glm::vec3(1.0f, 0.5f, 1.0f);
        //s->transform.Scale(scale);
        s->color = glm::vec3(1.0f, 0.0f, 0.0f);
        s->shadingType = Opaque;
        s->reflectivity = 0.85f;
        s->transmittance = 0.0f;
        s->IndexOfRefraction = 3.0f;
        objects.push_back((ObjectBase*)(s));
         
        // Add sphere
        Sphere* s2 = new Sphere();
        center = glm::vec3(0.f, -0.5f, -3.0f);
        s2->transform.Translate(center);
        scale = glm::vec3(0.5f, 0.25f, 0.5f);
        s2->transform.Scale(scale);
        s2->color = glm::vec3(0.0f, 1.0f, 0.0f);
        s2->shadingType = Opaque;
        s2->reflectivity = 0.0f;
        s2->transmittance = 0.0f;
        s2->IndexOfRefraction = 3.0;
        objects.push_back((ObjectBase*)(s2));

        Sphere* s3 = new Sphere();
        center = glm::vec3(2.0f, 0.5f, -5.0f);
        s3->transform.Translate(center);
        //scale = glm::vec3(1.0f, 0.5f, 1.0f);
        //s3->transform.Scale(scale);
        s3->color = glm::vec3(0.0f, 0.0f, 1.0f);
        s3->shadingType = Transparent;
        s3->reflectivity = 0.02f;
        s3->transmittance = 0.0f;
        s3->IndexOfRefraction = ior;
        objects.push_back((ObjectBase*)(s3));

        // Add Plane
        Plane* floor = new Plane();
        floor->planeNormal = glm::vec3(0.0f, 1.0f, 0.0f);
        floor->planePosition = glm::vec3(0.f, -1.0f, 0.0f);
        floor->useColorPattern = true;
        floor->color = glm::vec3(0.9f, 0.9f, 0.9f);
        floor->shadingType = Opaque;
        floor->reflectivity = 0.0f;
        floor->transmittance = 0.0f;
        floor->IndexOfRefraction = 1.0f;
        objects.push_back((ObjectBase*)(floor));

        // Add Plane
        Plane* back = new Plane();
        back->planeNormal = glm::vec3(0.0f, 0.f, 1.f);
        back->planeNormal = glm::normalize(back->planeNormal);
        back->planePosition = glm::vec3(0.f, 0.f, -11.f);
        back->useColorPattern = false;
        back->color = glm::vec3(0.9f, 0.9f, 0.9f);
        back->shadingType = Opaque;
        back->reflectivity = 0.02f;
        back->transmittance = 0.0f;
        back->IndexOfRefraction = 1.0f;
        objects.push_back((ObjectBase*)(back));

        // Add Plane
        Plane* left = new Plane();
        left->planeNormal = glm::vec3(1.0f, 0.0f, 0.0f);
        left->planePosition = glm::vec3(-4.0, 0.f, 0.0f);
        left->useColorPattern = false;
        left->color = glm::vec3(0.f, 0.f, 1.f);
        left->shadingType = Opaque;
        left->reflectivity = 0.02f;
        left->transmittance = 0.0f;
        left->IndexOfRefraction = 1.0f;
        objects.push_back((ObjectBase*)(left));

         // Add Plane
        Plane* right = new Plane();
        right->planeNormal = glm::vec3(-1.0f, 0.0f, 0.0f);
        right->planePosition = glm::vec3(4.0, 0.f, 0.0f);
        right->useColorPattern = false;
        right->color = glm::vec3(1.f, 0.0f, 0.0f);
        right->shadingType = Opaque;
        right->reflectivity = 0.02f;
        right->transmittance = 0.0f;
        right->IndexOfRefraction = 1.0f;
        objects.push_back((ObjectBase*)(right));
    }

    inline float DiffuseK(glm::vec3 lightdir, glm::vec3 pNorm)
    {
        return glm::clamp(glm::dot(lightdir, pNorm),-1.f, 0.f);
    }

    float SpecularK(glm::vec3 lightdir, glm::vec3 pHit, glm::vec3 pNorm, float a)
    {
        glm::vec3 viewDir = glm::normalize(rayTracer.camera.eyepos - pHit);
        glm::vec3 halfV = glm::normalize(viewDir + lightdir);
        float HdotN = glm::dot(halfV, pNorm);
        return pow(HdotN, a);
    }

    glm::vec3 BlinnPhong(const glm::vec3& color, const Light& light, const HitRecord& hitRec, const float specPower)
    {
        glm::vec3 lightdir = light.position  - hitRec.hitPoint;
        lightdir = glm::normalize(lightdir);
        float ambientK = 1.0f;
        float diffuseK = DiffuseK(lightdir, hitRec.hitNormal);
        float specularK = SpecularK(lightdir, hitRec.hitPoint, hitRec.hitNormal, specPower);
        return (ambientK * color) + (diffuseK * color) + (specularK * light.color);
    }

    string printVec(glm::vec3 v) const
    {
		stringstream ss;
		ss << "(" << v.x << "," << v.y << "," << v.z << ")";
		return ss.str();
    }

    public:
    std::vector<Light> lights;
    std::vector<ObjectBase*> objects;
    Renderer rayTracer;
    glm::vec3 offsets[4];

    private:
    bool initialized;
};