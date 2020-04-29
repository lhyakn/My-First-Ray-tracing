#ifndef SPHERE_H
#define SPHERE_H

#include "hittable.h"
#include "vec3.h"
#include "material.h"

class sphere: public hittable {
    public:
        sphere() {}
        sphere(vec3 cen, double r, shared_ptr<material> m)
            : center(cen), radius(r), mat_ptr(m),area(4 * pi * r * r),isemit(m->has_emition()) {};

        virtual bool hit(const ray& r, double tmin, double tmax, hit_record& rec) const;

        virtual aabb bounding_box() const {
            aabb output_box = aabb(
                center - vec3(radius, radius, radius),
                center + vec3(radius, radius, radius));
            return output_box;
        }

        virtual double getArea() const{
            return area;
        }

        virtual void Sample(hit_record& rec, double &pdf){
            float theta = 2.0 * pi * random_double(), phi = pi * random_double();
            vec3 dir(std::cos(phi), std::sin(phi)*std::cos(theta), std::sin(phi)*std::sin(theta));
            rec.p = center + radius * dir;
            rec.normal = dir;
            rec.mat_ptr = this->mat_ptr;
            pdf = 1.0d / area;
        }

        virtual double Isemit() const{
            return this->isemit;
        };

        virtual shared_ptr<material> get_mt() const{
            return mat_ptr;
        }

    public:
        vec3 center;
        double radius;
        shared_ptr<material> mat_ptr;
        double area;
        bool isemit;
};

#endif
