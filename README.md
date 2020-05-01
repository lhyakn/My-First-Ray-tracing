# My-First-Ray-tracing
参考了ray tracing one weekend的基于C++的光追渲染器，加上之前做的软光栅渲染器算是对图形学的一个总结。

# 实现功能
 1. **Path tracing**:Path tracing基本流程过了一遍。
 2. **Defocus blur**:模拟真实世界中相机景深的效果（大光圈怼脸）。
 3. **加速结构**:基于SAH构造的BVH加速结构。
 4. **多线程**:偷懒用了openMP实现的多线程，也算是实现了多线程叭。
 
# TODO LIST
 ~~1. 基于渲染方程的更加真实的path tracing，蒙特卡罗方法。~~
 1.1 微表面模型的实现(Microfacet Model)
 2. 更先进的光线追踪算法PPM之类的。

# 图片
  ![Balls](https://github.com/lhyakn/My-First-Ray-tracing/blob/master/bin/Release/Balls.png) 
  
  ![Cornell Box](https://github.com/lhyakn/My-First-Ray-tracing/blob/master/bin/Release/glass_box.png)
  
  ![Cornell Box PBR](https://github.com/lhyakn/My-First-Ray-tracing/blob/master/bin/Release/Cornell_box.png)
  采样数为1000的基于物理的全局光照渲染效果
  
# 参考
 1.[Ray Tracing in One Weekend](https://raytracing.github.io/) 
 
 2.[GAMES101-现代计算机图形学入门-闫令琪](https://www.bilibili.com/video/BV1X7411F744?p=22)
 
强烈安利GAMES101,闫老师讲的真的很好。
