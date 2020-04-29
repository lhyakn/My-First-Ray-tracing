#ifndef AABB_H
#define AABB_H

#include <limits>
#include "rtweekend.h"

class aabb {
    public:
        aabb() {
            double minNum = std::numeric_limits<double>::lowest();
            double maxNum = std::numeric_limits<double>::max();
            _max = vec3(minNum, minNum, minNum);
            _min = vec3(maxNum, maxNum, maxNum);
        }

        aabb(const vec3& a, const vec3& b) {
            _min = vec3(ffmin(a.x(), b.x()), ffmin(a.y(), b.y()), ffmin(a.z(), b.z()));
            _max = vec3(ffmax(a.x(), b.x()), ffmax(a.y(), b.y()), ffmax(a.z(), b.z()));
        }

        vec3 Diagonal() const { return _max - _min; }

        vec3 min() const {return _min; }
        vec3 max() const {return _max; }

        bool hit(const ray& r, double tmin, double tmax) const;

        vec3 _min;
        vec3 _max;

        int maxExtent() const
        {
            vec3 d = Diagonal();
            if (d.x() > d.y() && d.x() > d.z())
                return 0;
            else if (d.y() > d.z())
                return 1;
            else
                return 2;
        }

        double SurfaceArea() const
        {
            vec3 d = Diagonal();
            return 2 * (d.x() * d.y() + d.x() * d.z() + d.y() * d.z());
        }

        vec3 Centroid() { return 0.5 * _min + 0.5 * _max; }

        vec3 Offset(const vec3& p) const
        {
            vec3 o = p - _min;
            if (_max.x() > _min.x())
                o.e[0] /= _max.x() - _min.x();
            if (_max.y() > _min.y())
                o.e[1] /= _max.y() - _min.y();
            if (_max.z() > _min.z())
                o.e[2] /= _max.z() - _min.z();
            return o;
        }
};

inline bool aabb::hit(const ray& r, double tmin, double tmax) const {

    /*
    float t_min0,t_max0;
    float t_min1,t_max1;
    float t_min2,t_max2;
    vec3 unit_dir = unit_vector(r.direction());
    //x axis
    t_min0 = std::min((_min.x() - r.origin().x()) * (1/unit_dir.x()),(_max.x() - r.origin().x()) * (1/unit_dir.x()));
    t_max0 = std::max((_min.x() - r.origin().x()) * (1/unit_dir.x()),(_max.x() - r.origin().x()) * (1/unit_dir.x()));
    //y axis
    t_min1 = std::min((_min.y() - r.origin().y()) * (1/unit_dir.y()),(_max.y() - r.origin().y()) * (1/unit_dir.y()));
    t_max1 = std::max((_min.y() - r.origin().y()) * (1/unit_dir.y()),(_max.y() - r.origin().y()) * (1/unit_dir.y()));
    //z axis
    t_min2 = std::min((_min.z() - r.origin().z()) * (1/unit_dir.z()),(_max.z() - r.origin().z()) * (1/unit_dir.z()));
    t_max2 = std::max((_min.z() - r.origin().z()) * (1/unit_dir.z()),(_max.z() - r.origin().z()) * (1/unit_dir.z()));

    float t_enter = std::max(t_min0,std::max(t_min1,t_min2));
    float t_exit = std::min(t_max0,std::min(t_max1,t_max2));

    //std::cout << r.orig.x() << " " << r.orig.y() << " " << r.orig.z() << " " << std::endl;
    //std::cout << unit_dir.x() << " " << unit_dir.y() << " " << unit_dir.z() << " " << std::endl;

    //std::cout << t_enter << " " << t_exit <<  " " << (t_enter<t_exit&&t_exit>=0) << std::endl;
    //std::cout << std::endl;

    return (t_enter<=t_exit&&t_exit>=0);
    //else return false;
    */
     for (int a = 0; a < 3; a++) {
        auto invD = 1.0f / r.direction()[a];
        auto t0 = (min()[a] - r.origin()[a]) * invD;
        auto t1 = (max()[a] - r.origin()[a]) * invD;
        if (invD < 0.0f)
            std::swap(t0, t1);
        tmin = t0 > tmin ? t0 : tmin;
        tmax = t1 < tmax ? t1 : tmax;
        if (tmax < tmin)
            return false;
    }
    return true;
}

aabb surrounding_box(const aabb &box0,const aabb &box1);
aabb surrounding_box(const aabb &box0,const vec3 &p);

#endif // AABB_H
