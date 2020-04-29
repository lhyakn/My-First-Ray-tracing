#ifndef HITTABLE_LIST_H
#define HITTABLE_LIST_H

#include "hittable.h"
#include "bvh.h"
#include <memory>
#include <vector>

using std::shared_ptr;
using std::make_shared;

class hittable_list{
    public:
        hittable_list() {}
        hittable_list(shared_ptr<hittable> object) { add(object);}

        void clear() { objects.clear(); }
        void add(shared_ptr<hittable> object) { objects.push_back(object); area+=object->getArea();}

        bool hit(const ray& r, double tmin, double tmax, hit_record& rec) const;

        aabb bounding_box() const {
            aabb bounds;
            if(objects.size() == 0) return bounds;

            for (const auto& object : objects) {
                bounds = surrounding_box(bounds, object -> bounding_box());
            }

            return bounds;
        }

        double getArea() const{
            return this->area;
        }

        void buildBVH(){
            this->bvhl = new BVHAccel(objects, 1);
        }

        void sample_light(hit_record& rec,double &pdf) const;


    public:
        std::vector<shared_ptr<hittable>> objects;
        double area;
        BVHAccel *bvhl;
};



#endif
