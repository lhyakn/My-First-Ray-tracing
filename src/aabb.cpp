#include "ray.h"
#include "aabb.h"

aabb surrounding_box(const aabb &box0,const aabb &box1) {
    vec3 small(ffmin(box0.min().x(), box1.min().x()),
               ffmin(box0.min().y(), box1.min().y()),
               ffmin(box0.min().z(), box1.min().z()));
    vec3 big  (ffmax(box0.max().x(), box1.max().x()),
               ffmax(box0.max().y(), box1.max().y()),
               ffmax(box0.max().z(), box1.max().z()));
    return aabb(small,big);
}

aabb surrounding_box(const aabb &box0,const vec3 &p){
    vec3 small(ffmin(box0.min().x(), p.x()),
               ffmin(box0.min().y(), p.y()),
               ffmin(box0.min().z(), p.z()));
    vec3 big  (ffmax(box0.max().x(), p.x()),
               ffmax(box0.max().y(), p.y()),
               ffmax(box0.max().z(), p.z()));
    return aabb(small,big);
}

