#include "bvh.h"

#include <algorithm>
#include <cassert>
#include <iostream>

BVHAccel::BVHAccel(std::vector<shared_ptr<hittable>> p, int maxPrimsInNode)
    : maxPrimsInNode(std::min(255, maxPrimsInNode)),primitives(std::move(p))
{
    if (primitives.empty())
    {
        return;
    }

    root = recursiveBuild_SAH(primitives);
}


BVHBuildNode* BVHAccel::recursiveBuild_SAH(std::vector<shared_ptr<hittable>> objects)
{

    BVHBuildNode* node = new BVHBuildNode();
    aabb bounds;
    int objlen = objects.size();

    for (int i = 0; i < objlen; ++i)
        bounds = surrounding_box(bounds, objects[i]->bounding_box());
    if (objects.size() == 1) {
        // Create leaf _BVHBuildNode_
        node->bounds = objects[0]->bounding_box();
        node->object = objects[0];
        node->left = nullptr;
        node->right = nullptr;
        node->area = objects[0]->getArea();
        //std::cout << std::endl;
        //std::cout << node->bounds.min().x() <<  " " << node->bounds.min().y() <<  " " << node->bounds.min().z() << " " << std::endl;
        //std::cout << node->bounds.max().x() <<  " " << node->bounds.max().y() <<  " " << node->bounds.max().z() << " " << std::endl;
        return node;
    }
    else if (objlen == 2) {
        node->left = recursiveBuild_SAH(std::vector{objects[0]});
        node->right = recursiveBuild_SAH(std::vector{objects[1]});

        node->bounds = surrounding_box(node->left->bounds, node->right->bounds);
        node->area = node->left->area + node->right->area;
        return node;
    }
    else {
        aabb centroidBounds;
        for (int i = 0; i < objlen; ++i)
            centroidBounds = surrounding_box(centroidBounds, objects[i]->bounding_box().Centroid());
        int dim = centroidBounds.maxExtent();
        constexpr int nBuckets = 32;
        struct BucketInfo{
            int count = 0;
            aabb bounds;
        };
        float cost = std::numeric_limits<float>::max();
        int mincostSplitBucket = 0;
        BucketInfo buckets[nBuckets];
        //initialize
        for(int i = 0;i < objlen;++i){
            int b = nBuckets * centroidBounds.Offset(objects[i]->bounding_box().Centroid())[dim];
            if(b == nBuckets) b = nBuckets - 1;
            ++buckets[b].count;
            buckets[b].bounds = surrounding_box(buckets[b].bounds,objects[i]->bounding_box());
        }
        float tempcost = 0;
        for(int i = 0;i < nBuckets - 1;++i){
            aabb b0,b1;
            int count0=0,count1=0;
            for(int j=0;j<=i;++j){
                b0 = surrounding_box(b0,buckets[j].bounds);
                count0 += buckets[j].count;
            }
            for(int j=i+1;j<=nBuckets - 1;++j){
                b1 = surrounding_box(b1,buckets[j].bounds);
                count1 += buckets[j].count;
            }
            tempcost = 0.125f + (count0 * b0.SurfaceArea() + count1 * b1.SurfaceArea())/bounds.SurfaceArea();
            if(cost > tempcost){
                cost = tempcost;
                mincostSplitBucket = i;
            }
        }

        auto pmid = std::partition(objects.begin(),objects.end(),[=](const auto pi){
            int b = nBuckets * centroidBounds.Offset(pi->bounding_box().Centroid())[dim];
            if(b == nBuckets) b = nBuckets - 1;
            return b <= mincostSplitBucket;
        });
        auto beginning = objects.begin();
        auto middling = pmid;
        auto ending = objects.end();

        auto leftshapes = std::vector<shared_ptr<hittable>>(beginning, middling);
        auto rightshapes = std::vector<shared_ptr<hittable>>(middling, ending);

        assert(objects.size() == (leftshapes.size() + rightshapes.size()));

        node->left = recursiveBuild_SAH(leftshapes);
        node->right = recursiveBuild_SAH(rightshapes);
        //node->area = node->left->area + node->right->area;
        node->bounds = surrounding_box(node->left->bounds, node->right->bounds);
    }
    return node;


    /*
    //中点构造法
    BVHBuildNode* node = new BVHBuildNode();

    // Compute bounds of all primitives in BVH node
    aabb bounds;
    int objlen = objects.size();
    for (int i = 0; i < objlen; ++i)
        bounds = surrounding_box(bounds, objects[i]->bounding_box());
    if (objects.size() == 1) {
        // Create leaf _BVHBuildNode_
        node->bounds = objects[0]->bounding_box();
        node->object = objects[0];
        node->left = nullptr;
        node->right = nullptr;
        return node;
    }
    else if (objects.size() == 2) {
        node->left = recursiveBuild_SAH(std::vector{objects[0]});
        node->right = recursiveBuild_SAH(std::vector{objects[1]});

        node->bounds = surrounding_box(node->left->bounds, node->right->bounds);
        return node;
    }
    else {
        aabb centroidBounds;
        for (int i = 0; i < objlen; ++i)
            centroidBounds =
                surrounding_box(centroidBounds, objects[i]->bounding_box().Centroid());
        int dim = centroidBounds.maxExtent();
        float pmid = (centroidBounds.max()[dim] + centroidBounds.min()[dim])/2.0f;
        std::vector<shared_ptr<hittable>>::iterator it = std::partition(objects.begin(),objects.end(),[dim,pmid](shared_ptr<hittable> p1){
            return p1->bounding_box().Centroid()[dim] < pmid;
        });
        auto beginning = objects.begin();
        auto middling = it;
        auto ending = objects.end();
        auto leftshapes = std::vector<shared_ptr<hittable>>(beginning, middling);
        auto rightshapes = std::vector<shared_ptr<hittable>>(middling, ending);

        assert(objects.size() == (leftshapes.size() + rightshapes.size()));

        node->left = recursiveBuild_SAH(leftshapes);
        node->right = recursiveBuild_SAH(rightshapes);

        node->bounds = surrounding_box(node->left->bounds, node->right->bounds);
    }

    return node;
    */
}

bool BVHAccel::hit(const ray& r, double t_min, double t_max, hit_record& rec) const
{
    if (!root)
        return false;
    else
        return getIntersect(root,r,t_min,t_max,rec);
}

bool getIntersect(BVHBuildNode *node,const ray& r,double t_min, double t_max,hit_record& rec)
{

    if(node==nullptr)
        return false;
    else{
        if(node->bounds.hit(r,t_min,t_max)){
            if(node->left==nullptr&&node->right==nullptr){
                if(node->object->hit(r,t_min,t_max,rec))
                    return true;
                else
                    return false;
            }
            else{
                //hit_record temp_rec1;
                //temp_rec1.t = std::numeric_limits<double>::max();
                //hit_record temp_rec2;
                //temp_rec2.t = std::numeric_limits<double>::max();

                bool flag1 = getIntersect(node->left,r,t_min,t_max,rec);
                bool flag2 = getIntersect(node->right,r,t_min,flag1?rec.t:t_max,rec);


                return flag1||flag2;
            }
        }
        else
            return false;
    }
}

aabb BVHAccel::bounding_box() const
{
    return root ? root->bounds : aabb();
}

 void BVHAccel::getSample(BVHBuildNode* node, double p, hit_record& rec, double &pdf){
    if(node->left == nullptr || node->right == nullptr){
        node->object->Sample(rec, pdf);
        pdf *= node->area;
        return;
    }
    if(p < node->left->area) getSample(node->left, p, rec, pdf);
    else getSample(node->right, p - node->left->area, rec, pdf);
}

double BVHAccel::getArea() const{
    return root->area;
}

void BVHAccel::Sample(hit_record& rec,double &pdf){
    double p = std::sqrt(random_double()) * root->area;
    getSample(root, p, rec, pdf);
    pdf /= root->area;
}
