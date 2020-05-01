#ifndef MATERIAL_H
#define MATERIAL_H
#include "texture.h"

inline double schlick(double cosine, double ref_idx) {
    auto r0 = (1-ref_idx) / (1+ref_idx);
    r0 = r0*r0;
    return r0 + (1-r0)*pow((1 - cosine),5);
}

class material {
    public:
        virtual vec3 emitted(double u, double v, const vec3& p) const {
            return vec3(0,0,0);
        }

        virtual vec3 emitted() const {
            return vec3(0,0,0);
        }

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered
        ) const = 0;

        virtual vec3 eval(const vec3 &wi, const vec3 &wo, const vec3 &N) const = 0;

        virtual vec3 sample(const vec3 &wi, const vec3 &N) const = 0;

        virtual double pdf(const vec3 &wi,const vec3 &wo, const vec3 &N) const = 0;

        vec3 toWorld(const vec3 &a,const vec3 &N) const {
            vec3 B,C;
            if(std::fabs(N.x()) > std::fabs(N.y()))
            {
                double invLen = 1.0f / std::sqrt(N.x()*N.x()+N.z()*N.z());
                C = vec3(N.z() * invLen, 0.0f, -N.x() * invLen);
            }
            else
            {
                double invLen = 1.0f / std::sqrt(N.y()*N.y()+N.z()*N.z());
                C = vec3(0.0d,N.z() * invLen,-N.y() * invLen);
            }
            B = cross(C,N);
            return a.x() * B + a.y() * C + a.z() * N;
        }



    public:
        vec3 Kd;
};

class lambertian : public material {
    public:
        lambertian(shared_ptr<texture> a) : albedo(a) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered
        ) const {
            vec3 scatter_direction = rec.normal + random_unit_vector();
            scattered = ray(rec.p, scatter_direction);
            attenuation = albedo->value(rec.u, rec.v, rec.p);
            return true;
        }


        virtual vec3 eval(const vec3 &wi, const vec3 &wo, const vec3 &N) const {
            double cosalpha = dot(N, wo);
            if (cosalpha > 0.0d) {
                vec3 diffuse = Kd / pi;
                return diffuse;
            }
            else
                return vec3(0.0d,0.0d,0.0d);
        }

        virtual vec3 sample(const vec3 &wi, const vec3 &N) const{
            double x_1 = random_double(), x_2 = random_double();
            double z = std::fabs(1.0d - 2.0d * x_1);
            double r = std::sqrt(1.0d - z * z),phi = 2 * pi * x_2;
            vec3 localRay(r * std::cos(phi), r * std::sin(phi),z);
            return toWorld(localRay,N);
        }

        virtual double pdf(const vec3 &wi,const vec3 &wo, const vec3 &N) const{
            if(dot(wo,N) > 0.0d)
                return 0.5d / pi;
            else
                return 0.0f;
        }


    public:
        shared_ptr<texture> albedo;
};

class metal : public material {
    public:
        metal(shared_ptr<texture> a, double f) : albedo(a), fuzz(f < 1 ? f : 1) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered
        ) const {
            vec3 reflected = reflect(unit_vector(r_in.direction()), rec.normal);
            scattered = ray(rec.p, reflected + fuzz*random_in_unit_sphere());
            attenuation = albedo->value(rec.u, rec.v, rec.p);
            return (dot(scattered.direction(), rec.normal) > 0);
        }

        //待改
        virtual vec3 eval(const vec3 &wi, const vec3 &wo, const vec3 &N) const {
            double cosalpha = dot(N, wo);
            if (cosalpha > 0.0f) {
                vec3 diffuse = Kd / pi;
                return diffuse;
            }
            else
                return vec3(0.0f,0.0f,0.0f);
        }
        //待改
        virtual vec3 sample(const vec3 &wi, const vec3 &N) const{
            double x_1 = random_double(), x_2 = random_double();
            double z = std::fabs(1.0d - 2.0d * x_1);
            double r = std::sqrt(1.0d - z * z),phi = 2 * pi * x_2;
            vec3 localRay(r * std::cos(phi), r * std::sin(phi),z);
            return toWorld(localRay,N);
        }
        //待改
        virtual double pdf(const vec3 &wi,const vec3 &wo, const vec3 &N) const{
            if(dot(wo,N) > 0.0d)
                return 0.5d / pi;
            else
                return 0.0f;
        }

    public:
        shared_ptr<texture> albedo;
        double fuzz;
};

class dielectric : public material {
    public:
        dielectric(double ri) : ref_idx(ri) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered
        ) const {
            attenuation = vec3(1.0, 1.0, 1.0);
            double etai_over_etat;
            if (rec.front_face) {
                etai_over_etat = 1.0 / ref_idx;
            } else {
                etai_over_etat = ref_idx;
            }

            vec3 unit_direction = unit_vector(r_in.direction());
            double cos_theta = ffmin(dot(-unit_direction, rec.normal), 1.0);
            double sin_theta = sqrt(1.0 - cos_theta*cos_theta);
            if (etai_over_etat * sin_theta > 1.0 ) {
                vec3 reflected = reflect(unit_direction, rec.normal);
                scattered = ray(rec.p, reflected);
                return true;
            }

            double reflect_prob = schlick(cos_theta, etai_over_etat);
            if (random_double() < reflect_prob)
            {
                vec3 reflected = reflect(unit_direction, rec.normal);
                scattered = ray(rec.p, reflected);
                return true;
            }

            vec3 refracted = refract(unit_direction, rec.normal, etai_over_etat);
            scattered = ray(rec.p, refracted);
            return true;
        }

        //待改
        virtual vec3 eval(const vec3 &wi, const vec3 &wo, const vec3 &N) const {
            double cosalpha = dot(N, wo);
            if (cosalpha > 0.0f) {
                vec3 diffuse = Kd / pi;
                return diffuse;
            }
            else
                return vec3(0.0f,0.0f,0.0f);
        }
        //待改
        virtual vec3 sample(const vec3 &wi, const vec3 &N) const{
            double x_1 = random_double(), x_2 = random_double();
            double z = std::fabs(1.0d - 2.0d * x_1);
            double r = std::sqrt(1.0d - z * z),phi = 2 * pi * x_2;
            vec3 localRay(r * std::cos(phi), r * std::sin(phi),z);
            return toWorld(localRay,N);
        }
        //待改
        virtual double pdf(const vec3 &wi,const vec3 &wo, const vec3 &N) const{
            if(dot(wo,N) > 0.0d)
                return 0.5d / pi;
            else
                return 0.0f;
         }

        double ref_idx;
};

class diffuse_light : public material  {
    public:
        diffuse_light(shared_ptr<texture> a) : emit(a) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered
        ) const {
            return false;
        }

        virtual vec3 emitted(double u, double v, const vec3& p) const {
            return emit->value(u, v, p);
        }

        virtual vec3 emitted() const {
            return emit->value();
        }

        //待改
        virtual vec3 eval(const vec3 &wi, const vec3 &wo, const vec3 &N) const {
            double cosalpha = dot(N, wo);
            if (cosalpha > 0.0f) {
                vec3 diffuse = Kd / pi;
                return diffuse;
            }
            else
                return vec3(0.0f,0.0f,0.0f);
        }
        //待改
        virtual vec3 sample(const vec3 &wi, const vec3 &N) const{
            double x_1 = random_double(), x_2 = random_double();
            double z = std::fabs(1.0d - 2.0d * x_1);
            double r = std::sqrt(1.0d - z * z),phi = 2 * pi * x_2;
            vec3 localRay(r * std::cos(phi), r * std::sin(phi),z);
            return toWorld(localRay,N);
        }
        //待改
        virtual double pdf(const vec3 &wi,const vec3 &wo, const vec3 &N)const{
            if(dot(wo,N) > 0.0d)
                return 0.5d / pi;
            else
                return 0.0f;
         }

    public:
        shared_ptr<texture> emit;
};


#endif // MATERIAL_H
