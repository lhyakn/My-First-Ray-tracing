#include <ctime>
#include <iostream>
#include <math.h>


#include "rtweekend.h"
#include "hittable_list.h"
#include "sphere.h"
#include "camera.h"
#include "material.h"


#include "triangle.h"
#include "omp.h"



hittable_list random_scene() {
    hittable_list world;

    auto red = make_shared<lambertian>(make_shared<constant_texture>(vec3(0.65, 0.05, 0.05)));
    red->Kd = vec3(0.63d, 0.065d, 0.05d);
    auto white = make_shared<lambertian>(make_shared<constant_texture>(vec3(0.73, 0.73, 0.73)));
    white->Kd = vec3(0.725d, 0.71d, 0.68d);
    auto green = make_shared<lambertian>(make_shared<constant_texture>(vec3(0.12, 0.45, 0.15)));
    green->Kd = vec3(0.14d, 0.45d, 0.091d);
    auto light = make_shared<diffuse_light>(make_shared<constant_texture>(8.0f * vec3(0.747f+0.058f, 0.747f+0.258f, 0.747f) + 15.6f * vec3(0.740f+0.287f,0.740f+0.160f,0.740f) + 18.4f *vec3(0.737f+0.642f,0.737f+0.159f,0.737f)));
    light->Kd = vec3(0.65d,0.65d,0.65d);
    auto glass = make_shared<dielectric>(1.51630);

    //world.add(make_shared<flip_face>(make_shared<yz_rect>(0, 555, 0, 555, 555, green)));
    //world.add(make_shared<yz_rect>(0, 555, 0, 555, 0, red));
    //world.add(make_shared<xz_rect>(213, 343, 227, 332, 554, light));
    //world.add(make_shared<flip_face>(make_shared<xz_rect>(0, 555, 0, 555, 555, white)));
    //world.add(make_shared<xz_rect>(0, 555, 0, 555, 0, white));
    //world.add(make_shared<flip_face>(make_shared<xy_rect>(0, 555, 0, 555, 555, white)));

    //MeshTriangle floor("E:/ray tracer/Ray tracer/models/cornellbox/floor.obj", white);
    //MeshTriangle shortbox("E:/ray tracer/Ray tracer/models/cornellbox/shortbox.obj", white);
    //MeshTriangle tallbox("E:/ray tracer/Ray tracer/models/cornellbox/tallbox.obj", white);
    //MeshTriangle left("E:/ray tracer/Ray tracer/models/cornellbox/left.obj", red);
    //MeshTriangle right("E:/ray tracer/Ray tracer/models/cornellbox/right.obj", green);
    //MeshTriangle light_("E:/ray tracer/Ray tracer/models/cornellbox/light.obj", light);

    world.add(make_shared<MeshTriangle>("E:/ray tracer/Ray tracer/models/cornellbox/floor.obj",white));
    world.add(make_shared<MeshTriangle>("E:/ray tracer/Ray tracer/models/cornellbox/shortbox.obj", white));
    world.add(make_shared<MeshTriangle>("E:/ray tracer/Ray tracer/models/cornellbox/tallbox.obj", white));
    world.add(make_shared<MeshTriangle>("E:/ray tracer/Ray tracer/models/cornellbox/left.obj", red));
    world.add(make_shared<MeshTriangle>("E:/ray tracer/Ray tracer/models/cornellbox/right.obj", green));
    world.add(make_shared<MeshTriangle>("E:/ray tracer/Ray tracer/models/cornellbox/light.obj", light));

    return world;
}

vec3 ray_color(const ray& r, const vec3& background, const hittable_list& world, int depth, double RussianRoulette) {
    /*
    hit_record rec;
    hit_record pos;
    double pdf;
    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return vec3(0,0,0);

    // If the ray hits nothing, return the background color.
    if (!world.hit(r, 0.0001, infinity, rec))
        return background;

    ray scattered;
    vec3 attenuation;
    vec3 emitted = rec.mat_ptr->emitted(rec.u, rec.v, rec.p);
    // 如果打到了光源
    if (!rec.mat_ptr->scatter(r, rec, attenuation, scattered))
        return emitted;

    // 对光源采样


    world.sample_light(pos,pdf);
    vec3 position = pos.p;
    vec3 wi = unit_vector(r.direction());
    vec3 wo = unit_vector(position - rec.p);
    vec3 NN = pos.normal;
    ray temp(rec.p,pos.p-rec.p);
    emitted = pos.mat_ptr->emitted();

    if(!world.hit(temp, 0.0001, infinity, pos)){
            if(pos.p-)
        //emitted = emitted * rec.mat_ptr->eval(wi,wo,rec.normal) * dot(wo,rec.normal) * dot(wo,NN) / (position-rec.p).length_squared() / pdf;
    }
    else emitted = emitted - emitted;


    return emitted + attenuation * ray_color(scattered, background, world, depth-1);
    */
    hit_record intersect;
    hit_record lightpos;
    hit_record blockcheck;
    double pdf_light;
    vec3 L_dir(0.0d,0.0d,0.0d);
    vec3 L_indir(0.0d,0.0d,0.0d);
    vec3 hitcolor = background;


    if (depth <= 0)
        return vec3(0,0,0);

    //如果没有击中
    if (!world.hit(r, 0.0001, infinity, intersect))
        return hitcolor;

    vec3 emitted = intersect.mat_ptr->emitted();

    //如果击中了光源
    if (emitted.length_squared() > EPSILON)
        return emitted;


    //对光源进行采样
    world.sample_light(lightpos,pdf_light);
    vec3 wo = unit_vector(lightpos.p - intersect.p);
    vec3 wi = unit_vector(r.direction());
    vec3 N = unit_vector(intersect.normal);
    ray p2light(intersect.p,lightpos.p - intersect.p);

    if(world.hit(p2light,0.0001,infinity,blockcheck))
    {
        //如果中间没有被阻挡
        if(std::fabs((lightpos.p - intersect.p).length() - (blockcheck.p - intersect.p).length()) < EPSILON)
        {
            vec3 NN = unit_vector(blockcheck.normal);
            vec3 emit = lightpos.mat_ptr->emitted();
            L_dir = emit * intersect.mat_ptr->eval(wi,wo,N) * dot(wo,N) * dot(-wo,NN) / (lightpos.p-intersect.p).length_squared() / pdf_light;
        }
    }

    //计算间接光照//俄罗斯轮盘赌


    double P_RR = random_double();

    if( P_RR < RussianRoulette )
    {
        wo = unit_vector(intersect.mat_ptr->sample(wi,N));
        ray bounce(intersect.p,wo);
        hit_record tracebounce;
        if(world.hit(bounce,0.001,infinity,tracebounce))
        {
            if(tracebounce.mat_ptr->emitted().length() < EPSILON)
            {
                L_indir = ray_color(bounce,background,world,depth-1,RussianRoulette) * intersect.mat_ptr->eval(wi,wo,N) * dot(wo,N) / intersect.mat_ptr->pdf(wi,wo,N) / RussianRoulette;
            }
        }
    }

    hitcolor = L_dir + L_indir;
    return hitcolor;
}



int main() {
    const int image_width = 1024;
    const int image_height = 1024;
    const int samples_per_pixel = 1000;
    const int max_depth = 50;
    const auto aspect_ratio = double(image_width) / image_height;
    const vec3 background(0,0,0);
    const double RussianRoulette = 0.8d;

    auto world = random_scene();

    world.buildBVH();

    vec3 lookfrom(278, 273, -800);
    vec3 lookat(278,273,0);
    vec3 vup(0,1,0);
    auto dist_to_focus = 10.0;
    auto aperture = 0.0;
    auto vfov = 40.0;

    camera cam(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, dist_to_focus);
    std::cout << "P3\n" << image_width << " " << image_height << "\n255\n";

    std::vector<std::vector<vec3> > framebuffer(image_height,std::vector<vec3>(image_width));

    clock_t start_time=clock();
    for (int j = image_height-1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        #pragma omp parallel for
        for (int i = 0; i < image_width; ++i) {
            vec3 color(0, 0, 0);
            for (int s = 0; s < samples_per_pixel; ++s) {
                auto u = (i + random_double()) / image_width;
                auto v = (j + random_double()) / image_height;
                ray r = cam.get_ray(u, v);
                color += ray_color(r, background, world, max_depth, RussianRoulette);
            }
            framebuffer[j][i] = std::move(color);
        }
    }

    for (int j = image_height-1; j >= 0; --j){
        for(int i = 0; i < image_width; ++i){
            framebuffer[j][i].write_color(std::cout, samples_per_pixel);
        }
    }

    clock_t end_time=clock();

    std::cerr << "\nThe run time is: " <<(double)(end_time - start_time) / CLOCKS_PER_SEC << "s" ;
    std::cerr << "\nDone.\n";
}
