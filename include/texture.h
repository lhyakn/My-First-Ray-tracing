#ifndef TEXTURE_H
#define TEXTURE_H


#include "rtweekend.h"
#include "vec3.h"

class texture {
    public:
        virtual vec3 value(double u, double v, const vec3& p) const = 0;
        virtual vec3 value() const = 0;
};

class constant_texture : public texture {
    public:
        constant_texture() {}
        constant_texture(vec3 c) : color(c) {}

        virtual vec3 value(double u, double v, const vec3& p) const {
            return color;
        }
        virtual vec3 value() const {
            return color;
        }

    public:
        vec3 color;
};

#endif // TEXTURE_H
