# ray_tracer
A **learning** repo for https://raytracing.github.io/

ALL HONOR BELONGS TO Peter Shirley.

环境: Intel E2630 v4

# 目标
- [x] 完成三本ray tracer
- [x] cpu并行/~~移植到cuda上~~

# 其他
**一点经验**: 

- 并行很重要,越早实现越好。图形学里面debugger基本是残疾状态，除了诊断法向量之类的有点用。没有并行调试一次等半天。尤其是涉及到采样上的Bug, 采样数量不上来看不到Bug源头，没有实现并行效率会很低。

- hard-code参数不可取。应该参数和代码分离。这里使用了ini作为配置文件，表达力比较弱, 不能嵌套结构。如果从头再来,我会使用有嵌套树状结构的配置文件(json, toml等，可以在配置文件里面定义渲染场景，场景和渲染部分解耦。hard-code参数会导致频繁重新编译代码。

# 渲染图片
## cornell box(蒙特卡罗积分)
**image**目录有原图
- 金属、玻璃球(1200x1200p, 2048samples,时长941秒)

![cbox1](image/9.jpg =600x600)
- 传统cornell box

![cbox2](image/8.jpg =600x600)

## 其他有意思的结果/场景

- [x] ray tracing in one weekend(cpu并行,1600p x 1600p x 100 samples,渲染时长: 290s)
![one weekend](image/1.jpg =600x600)

- [x] Bounding Volume Hierarchies优化 (CPU并行，1600p x 1600p x 300 samples,渲染时长: 80.1s)

![BVH](image/2.jpg =600x600)

- [x] perlin noise(柏林噪音)

![perlin1](image/3.jpg =600x600)

![perlin2](image/4.jpg =600x600)

- [x] sphere uv mapping

![earth](image/5.jpg =600x600)

- [x] cornell box(没有用Monte-Carlo,随机均匀采样，噪点有点严重)

![box1](image/6.jpg =600x600)

![box2](image/7.jpg =600x600)
