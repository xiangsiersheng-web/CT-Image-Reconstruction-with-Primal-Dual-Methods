# 原始对偶方法求解CT图像重建问题

## 前言

这是本科时做的一个毕业设计，想要使用原始对偶方法求解一个欠定方程组，参考了此前的PDHG算法，后面又使用了变步长的思想。

其中原始对偶方法PDHG以及其变体CPPD、HYPD（隐式梯度下降）是参考了这三篇论文：

- 1.PDHG: ZHU M, CHAN T. An efficient primal-dual hybrid gradient algorithm for total variation image restoration[J]. Ucla Cam Report, 2008, 34: 8-34.

- 2.CPPD: CHAMBOLLE A, POCK T. A first-order primal-dual algorithm for convex problems with applications to imaging[J]. Journal of Mathematical Imaging and Vision, 2011, 40: 120-145

- 3.HYPD: HE B, YUAN X. Convergence analysis of primal-dual algorithms for a saddle-point problem: from contraction perspective[J]. SIAM Journal on Imaging Sciences, 2012, 5(1): 119-149.

变步长（自适应步长）思想来源于这篇文章：

- Goldstein T , Li M , Yuan X , et al. Adaptive Primal-Dual Hybrid Gradient Methods for Saddle-Point Problems: , 10.48550/arXiv.1305.0546[P]. 2013.

## 文件目录

- systemMatrix：生成2D图片的系统投影矩阵

- art：两种代数重建技术ART/SART的实现
- PD-rof：原始对偶方法HYPD求解rof问题的实现
- PD-TVCT：原始对偶方法HYPD/CPPD求解TVCT问题的实现（含变步长方法）
- 源图：绘制的图片
- compare_cppd&sart.mlx：将HYPD算法与CPPD、SART算法比较，凸显HYPD算法的优势，以及变步长的加速收敛
- parameter_sensitivity.mlx：参数灵敏度分析