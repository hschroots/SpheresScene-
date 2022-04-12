#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include <exception>
#include <vector>
#include <numeric>

#include "Renderer.h"
#include "Camera.h"
#include "ObjectBase.h"
#include "Sphere.h"
#include "Plane.h"
#include "Light.h"
#include "Ray.h"

#define MAX_TRACE_DEPTH 4

glm::vec3 clearColor(0.2f, 0.8f, 0.9f);

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
                float i = col * dx;
                float j = row * dy;

                rayTracer.computeRay(i, j, ray);

                glm::vec3 colorVal = glm::vec3(0); 
                for(int r=0; r < 4; r++)
                {
                    if( col == 200 && row == 238)
                    {
                        debug = true;
                        std::cout << "(row,col, r) = (" << row << "," << col << "," << r << ")" << std::endl;
                        std::cout << "Ray: " << ray.direction.x << ", " << ray.direction.y << ", " << ray.direction.z << std::endl;
                    }
                    ray.direction = glm::normalize(ray.direction + offsets[r]);
                    colorVal += 0.25f * traceRay(ray, 0, debug);
                }
                colorVal = glm::clamp(colorVal, 0.0f, 1.0f);

                const float b = colorVal.b;
                const float g = colorVal.g;
                const float r = colorVal.r;

                rayTracer.framebuffer.SetPixel(Pixel(b,g,r), col, row);
            }
        std::stringstream filename;
        filename << "out/images/output_ior_" << objects[0]->IndexOfRefraction << "_" << rayTracer.width << "x" << rayTracer.height << ".bmp";
        rayTracer.framebuffer.Export(filename.str().c_str());
        }
    }
    void fresnel(glm::vec3& incidentRay, glm::vec3& norm, float& ior, float& kr, bool debug = false)
    {
        glm::vec3 I = glm::normalize(incidentRay);
        glm::vec3 N = glm::normalize(norm);

        float cosi = glm::clamp(-1.0f, 1.0f, glm::dot(N, I));
        float etai = 1.0f;
        float etat = ior;
        
        if(cosi > 0.0f)
        {
            std::swap(etai, etat);
        }
        float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi)); 
        // Total internal reflection
        if (sint >= 1.0f) { 
            kr = 1.0f; 
        } 
        else { 
            float cost = sqrtf(std::max(0.f, 1.f - sint * sint)); 
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
            return glm::vec3(0.0f);
        }
        else
        {
            return eta * I + ( eta * cosi - std::sqrt(k) ) * N;
        }
    }

    glm::vec3 reflect(glm::vec3 incidentRay, glm::vec3 norm, bool debug = false)
    {
        glm::vec3 I = glm::normalize(incidentRay);
        glm::vec3 N = glm::normalize(norm);

        return 2.f * glm::dot(I, N) * N - I;
    }
    float traceShadowRay(Ray& ray, bool debug)
    {
        for(int idx = 0; idx < objects.size(); idx++)
        {
            HitRecord rec;
            if(rayTracer.intersect(ray, objects[idx], rec, debug))
            {
                return 1.0f;
            }
        }
        return 0.0f;
    }

    glm::vec3 traceRay(Ray& ray, uint32_t traceDepth, bool debug = false)
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
            if(rayTracer.intersect(ray, objects[idx], rec, debug))
            {
                float distance = glm::distance(ray.origin, rec.hitPoint);
                if(debug)
                {
                    std::cout << "Intersected with object idx =  " << idx << std::endl;
                    std::cout << "Object is at a distance of  " << distance << std::endl;
                    std::cout << "Is that less than the current distance of " << minDist << std::endl;
                }
                if(distance < minDist)
                {
                    obj = objects[idx];
                    hitRec = rec;
                    if (glm::length(hitRec.hitNormal) < 0.99999f)
                    {
                        hitRec.print();
                        std::cout << "Object Idx = " << idx << std::endl;
                    }
                    minDist = distance;
                }
            }
        }

        // If the min distance is still some large number, then the ray didn't hit anything
        // Set the color to the background color
        if(minDist >= 9999999999.99f)
        {
            return clearColor;
        }
        else // we hit an object and it is the closest object
        { 
            //Use only the first light in list
            Light light = lights[0];
            // trace shadow rays
            glm::vec3 lightdir = glm::normalize(light.position - hitRec.hitPoint);
            //float shadowK = 0.0f;
            //int numShadowRays = 4;
            //float invShadowRays = 1.f/(float)numShadowRays;
            //
            //for( int i = 0; i < numShadowRays; i++)
            //{
            //    Ray shadowRay;
            //    shadowRay.origin = worldPositionHit + (0.0001f * normalHit);
            //    shadowRay.direction = glm::normalize(lightdir + offsets[i % 4]);
            //    shadowK += invShadowRays * traceShadowRay(shadowRay, debug);
            //}

            switch(obj ? obj->shadingType : -1)
            {
                case BlinnPhong: // If the object is opaque shade it according to the blinn-phong model
                {
                    float ambientK = 0.4f;
                    float diffuseK = DiffuseK(lightdir, hitRec.hitNormal);
                    float specularK = SpecularK(lightdir, hitRec.hitPoint, hitRec.hitNormal, 100.f);

                    glm::vec3 specularColor(1.0f, 1.0f, 1.0f);
                    color = (ambientK * obj->color) + (diffuseK * obj->color) + (specularK * obj->color);
                }
                break;
                case ReflectAndRefract: // If the object is transparent trace reflection and refraction rays
                {
                    glm::vec3 refractColor(0.f);
                    glm::vec3 reflectColor(0.f);
                    float kr = 0.0f;
                    // if fresnel number is 1, then we have total internal reflection
                    fresnel(ray.direction, hitRec.hitNormal, obj->IndexOfRefraction, kr, debug);
                    bool outside = glm::dot(ray.direction, hitRec.hitNormal) < 0.f;
                    glm::vec3 bias = 0.0001f * hitRec.hitNormal;

                    if( kr < 1.0f) // trace refraction
                    {
                        Ray refractRay;
                        refractRay.origin = hitRec.front_face ? hitRec.hitPoint - bias : hitRec.hitPoint + bias;
                        refractRay.direction = glm::refract(ray.direction, hitRec.hitNormal, obj->IndexOfRefraction);
                        //trace refracted ray
                        refractColor = traceRay(refractRay, traceDepth +1, debug);
                    }
                    //trace refelction
                    Ray reflectRay;
                    reflectRay.origin = hitRec.front_face ? hitRec.hitPoint + bias : hitRec.hitPoint - bias;
                    reflectRay.direction = glm::reflect(ray.direction, hitRec.hitNormal);
                    reflectColor = traceRay(reflectRay, traceDepth + 1, debug);

                    color += reflectColor * kr + refractColor * (1 - kr);

                    //if(traceDepth == 0)
                    //{
                    //    float ambientK = 1.0f;
                    //    float diffuseK = DiffuseK(lightdir, normalHit);
                    //    float specularK = SpecularK(lightdir, worldPositionHit, normalHit, 100.f);

                    //    color += (diffuseK * color) + (specularK * color);
                    //    color = glm::clamp(glm::vec3(0.0f), glm::vec3(1.0f), color);
                    //}
                }
                break;
                default:
                break;
            }
            return color = hitRec.hitNormal;// *(1.f - (0.4f * shadowK));
        }
    }

    private:
    void setupLights()
    {
        Light light;
        light.position = glm::vec3(-2.5f, 3.5f, -3.0f);
        light.color = glm::vec3(1.0f, 1.0f, 1.0f);
        lights.push_back(light);
    }

    void setupGeometry(float ior = 1.1f)
    {
        // Add sphere
        Sphere* s = new Sphere();
        //glm::vec3 scale = glm::vec3(1.0f, 0.5f, 1.0f);
        //s->transform.Scale(scale);
        glm::vec3 center = glm::vec3(0.6, 0.5f, 4.0f);
        s->transform.Translate(center);
        s->color = glm::vec3(0.7f, 0.5f, 0.0f);
        s->shadingType = ReflectAndRefract;
        s->IndexOfRefraction = ior;
        objects.push_back((ObjectBase*)(s));

        // Add sphere
        Sphere* s2 = new Sphere();
        center = glm::vec3(-1.2f, 0.0f, 5.0f);
        s2->transform.Translate(center);
        s2->color = glm::vec3(0.7f, 0.5f, 0.0f);
        s2->shadingType = BlinnPhong;
        s2->IndexOfRefraction = 1.0f;
        objects.push_back((ObjectBase*)(s2));

        // Add Plane
        Plane* p = new Plane();
        p->planeNormal = glm::vec3(0.0f, 1.0f, 0.0f);
        p->planePosition = glm::vec3(0.f, -2.0f, 1.0f);
        p->color = glm::vec3(0.3f, 0.1f, 0.9f);
        p->shadingType = BlinnPhong;
        p->IndexOfRefraction = 1.0f;
        objects.push_back((ObjectBase*)(p));
    }

    inline float DiffuseK(glm::vec3 lightdir, glm::vec3 pNorm)
    {
        return std::max(glm::dot(lightdir, pNorm),0.f);
    }

    float SpecularK(glm::vec3 lightdir, glm::vec3 pHit, glm::vec3 pNorm, float a)
    {
        glm::vec3 viewDir = glm::normalize(rayTracer.camera.eyepos - pHit);
        glm::vec3 halfV = glm::normalize(viewDir + lightdir);
        float HdotN = glm::dot(halfV, pNorm);
        return pow(HdotN, a);
    }

    public:
    std::vector<Light> lights;
    std::vector<ObjectBase*> objects;
    Renderer rayTracer;
    glm::vec3 offsets[4];

    private:
    bool initialized;
};