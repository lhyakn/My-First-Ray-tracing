#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "hittable.h"
#include "vec3.h"
#include "OBJ_Loader.hpp"
#include "bvh.h"

#include <array>
#include <assert.h>
#include <algorithm>

class triangle : public hittable
{
    public:
        triangle();

        triangle(vec3 _v0, vec3 _v1, vec3 _v2, shared_ptr<material> _m)
        : v0(_v0), v1(_v1), v2(_v2), mat_ptr(_m)
        {
            e1 = v1 - v0;
            e2 = v2 - v0;
            normal = unit_vector(cross(e1,e2));
            area = cross(e1,e2).length()*0.5f;
        }

        virtual bool hit(const ray& r, double tmin, double tmax, hit_record& rec) const{
            vec3 unit_direction = unit_vector(r.direction());
            if(dot(unit_direction,normal) > 0)
                return false;
            double u,v;
            vec3 pvec = cross(unit_direction,e2);
            double det = dot(e1,pvec);
            if(std::fabs(det) < 0.00001)
                return false;
            double det_inv = 1. / det;
            vec3 tvec = r.origin() - v0;
            u = dot(tvec,pvec) * det_inv;
            if(u < 0 || u > 1)
                return false;
            vec3 qvec = cross(tvec,e1);
            v = dot(unit_direction,qvec) * det_inv;
            if(v < 0 || u + v > 1)
                return false;



            vec3 E1 = v1 - v0;
            vec3 E2 = v2 - v0;
            vec3 S  = r.origin() - v0;
            vec3 S1 = cross(unit_direction,E2);
            vec3 S2 = cross(S,E1);
            vec3 temp(dot(S2,E2),dot(S1,S),dot(S2,unit_direction));
            vec3 tuv = 1/dot(S1,E1) * temp;
            if(tuv.e[0]>tmin&&tuv.e[0]<tmax&&tuv.e[1]>0&&tuv.e[2]>0&&(1-tuv.e[1]-tuv.e[2])>0){
                rec.p = vec3(r.origin() + unit_direction * tuv.e[0]);
                rec.normal = normal;
                rec.mat_ptr = mat_ptr;
                rec.set_face_normal(r, normal);
                rec.t = tuv.e[0];
                return true;
            }
            return false;

        }

        virtual aabb bounding_box() const {
            aabb output_box =  surrounding_box(aabb(v0, v1),v2);
            return output_box;
        }

        virtual double getArea() const{
            return this -> area;
        }

        virtual void Sample(hit_record& rec,double &pdf){
            double x = std::sqrt(random_double()), y = random_double();
            rec.p = v0 * (1.0f - x) + v1 * (x * (1.0f - y)) + v2 * (x * y);
            rec.normal = this->normal;
            rec.mat_ptr = this->mat_ptr;
            pdf = 1.0f / area;
        }

        virtual double Isemit() const{
            return this->isemit;
        };

        virtual shared_ptr<material> get_mt() const{
            return mat_ptr;
        }

    public:
        vec3 v0,v1,v2;  //vertices A,B,C
        vec3 e1,e2;     //weo edges v1-v0,v2-v0;
        vec3 t0, t1, t2; // texture coords
        vec3 normal;
        double area;
        bool isemit;
        shared_ptr<material> mat_ptr;
};


class MeshTriangle : public hittable
{
    public:
        MeshTriangle(const std::string& filename, shared_ptr<material> mt)
        {
            objl::Loader loader;
            loader.LoadFile(filename);
            area = 0;
            mat_ptr = mt;
            assert(loader.LoadedMeshes.size() == 1);
            auto mesh = loader.LoadedMeshes[0];

            vec3 min_vert = vec3{std::numeric_limits<double>::infinity(),
                                         std::numeric_limits<double>::infinity(),
                                         std::numeric_limits<double>::infinity()};
            vec3 max_vert = vec3{-std::numeric_limits<double>::infinity(),
                                         -std::numeric_limits<double>::infinity(),
                                         -std::numeric_limits<double>::infinity()};
            int lenmesh = mesh.Vertices.size();
            for (int i = 0; i < lenmesh; i += 3) {
                std::array<vec3, 3> face_vertices;

                for (int j = 0; j < 3; j++) {
                    auto vert = vec3(mesh.Vertices[i + j].Position.X,
                                         mesh.Vertices[i + j].Position.Y,
                                         mesh.Vertices[i + j].Position.Z);
                    face_vertices[j] = vert;

                    min_vert = vec3(std::min(min_vert.e[0], vert.e[0]),
                                    std::min(min_vert.e[1], vert.e[1]),
                                    std::min(min_vert.e[2], vert.e[2]));

                    max_vert = vec3(std::max(max_vert.e[0], vert.e[0]),
                                    std::max(max_vert.e[1], vert.e[1]),
                                    std::max(max_vert.e[2], vert.e[2]));
                }

                triangles.emplace_back(face_vertices[0], face_vertices[1],
                                       face_vertices[2], mt);
            }

            bbox = aabb(min_vert, max_vert);

            std::vector<std::shared_ptr<hittable> > ptrs;

            for (auto& tri : triangles){
                ptrs.emplace_back(make_shared<triangle>(tri.v0,tri.v1,tri.v2,tri.mat_ptr));
                area += tri.getArea();
            }

            bvh = new BVHAccel(ptrs,1);
        }

        virtual bool hit(const ray& r, double tmin, double tmax, hit_record& rec) const{
            /*
            hit_record temp_rec;
            bool hit_anything = false;
            auto closest_so_far = tmax;

            for (auto& tri : triangles){
                if(tri.hit(r,tmin,closest_so_far,temp_rec))
                {
                    hit_anything = true;
                    closest_so_far = rec.t;
                    rec = temp_rec;
                }
            }

            return hit_anything;*/
            return this->bvh->hit(r,tmin,tmax,rec);
        }

        virtual aabb bounding_box() const {
            return bbox;
        }

        virtual double getArea() const {
            return this -> area;
        }

        virtual void Sample(hit_record& rec,double &pdf){
            this->bvh->Sample(rec, pdf);
            rec.mat_ptr = mat_ptr;
        }

        virtual double Isemit() const{
            return this->isemit;
        };

        virtual shared_ptr<material> get_mt() const{
            return mat_ptr;
        }
    public:
        std::shared_ptr<vec3[]> vertices;
        uint32_t numTriangles;
        std::shared_ptr<uint32_t[]> vertexIndex;
        std::vector<triangle> triangles;
        double area;
        bool isemit;
        shared_ptr<material> mat_ptr;
        aabb bbox;
        BVHAccel* bvh;
};

#endif // TRIANGLE_H
