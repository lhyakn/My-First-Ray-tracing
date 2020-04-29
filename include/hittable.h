#ifndef HITTABLE_H
#define HITTABLE_H


#include "ray.h"
#include "aabb.h"

class material;

struct hit_record {
    vec3 p;
    vec3 normal;
    shared_ptr<material> mat_ptr;
    double t;
    double u;
    double v;
    bool front_face;

    inline void set_face_normal(const ray& r, const vec3& outward_normal) {
        front_face = dot(r.direction(), outward_normal) < 0;
        normal = front_face ? outward_normal :-outward_normal;
    }
};

class hittable {
    public:
        virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const = 0;
        virtual aabb bounding_box() const = 0;
        virtual double getArea() const = 0;
        virtual void Sample(hit_record& rec,double &pdf) = 0;
        virtual double Isemit() const{
            return this->isemit;
        };
        virtual shared_ptr<material> get_mt() const{
            return mat_ptr;
        }
    public:
        double area;
        bool isemit ;
        shared_ptr<material> mat_ptr;
};

/*
class flip_face : public hittable {
    public:
        flip_face(shared_ptr<hittable> p) : ptr(p) {}

        virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
            if (!ptr->hit(r, t_min, t_max, rec))
                return false;

            rec.front_face = !rec.front_face;
            return true;
        }


    public:
        shared_ptr<hittable> ptr;
};
*/
#endif
