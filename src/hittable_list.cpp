#include "hittable_list.h"
#include "material.h"
bool hittable_list::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    /*
    hit_record temp_rec;

    bool hit_anything = false;
    auto closest_so_far = t_max;


    for (const auto& object : objects) {
        if (object->hit(r, t_min, closest_so_far, temp_rec)) {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            rec = temp_rec;
        }
    }
    */

    bool flag = this->bvhl->hit(r,t_min,t_max,rec);

    return flag;
}

void hittable_list::sample_light(hit_record& rec,double &pdf) const{
    double total_area = 0;
    for (uint32_t k = 0; k < objects.size(); ++k){
        if(objects[k]->Isemit()){
            total_area += objects[k]->getArea();
        }
    }
    //std::cout << total_area << std::endl;

    double p = random_double() * total_area;
    total_area = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->Isemit()){
            total_area += objects[k]->getArea();
            if (p <= total_area){
                objects[k]->Sample(rec, pdf);
                break;
            }
        }
    }
}
