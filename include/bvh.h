#ifndef BVH_H
#define BVH_H

#include <memory>
#include <vector>
#include "hittable.h"

using std::shared_ptr;
using std::make_shared;

struct BVHBuildNode;
// BVHAccel Forward Declarations
struct BVHPrimitiveInfo;

class BVHAccel : public hittable {

public:

    // BVHAccel Public Methods
    BVHAccel(std::vector<shared_ptr<hittable>> p, int maxPrimsInNode = 1);
    ~BVHAccel();

    virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const;

    virtual aabb bounding_box() const;

    virtual double getArea() const;

    virtual void Sample(hit_record& rec,double &pdf);

    void getSample(BVHBuildNode* node, double p, hit_record& rec, double &pdf);
    // SAH
    BVHBuildNode* recursiveBuild_SAH(std::vector<shared_ptr<hittable>>objects);

public:
    BVHBuildNode* root;
    // BVHAccel Private Data
    const int maxPrimsInNode;

    std::vector<shared_ptr<hittable>> primitives;
};

struct BVHBuildNode {
    aabb bounds;
    double area;
    BVHBuildNode *left;
    BVHBuildNode *right;
    shared_ptr<hittable> object;

public:
    int splitAxis=0, firstPrimOffset=0, nPrimitives=0;
    // BVHBuildNode Public Methods
    BVHBuildNode(){
        bounds = aabb();
        left = nullptr;right = nullptr;
        object = nullptr;
    }
};
bool getIntersect(BVHBuildNode *node,const ray& r,double t_min, double t_max,hit_record& rec);
#endif // BVH_H
