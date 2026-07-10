# 中心网 (CenterNetworkECG)

从一组网络中提取最具代表性的**中心网络**。

**GitHub:** [Yu-Wang-0923/CenterNetworkECG](https://github.com/Yu-Wang-0923/CenterNetworkECG)

---

## EMD 是什么

EMD (Edgewise Multiplex Depth) 是一篇 2025 年发表在 JCGS 上的统计方法论文（Zhang, Deng & Dai），提出了一种针对网络型数据（network-valued data）的统计深度（statistical depth）概念。

核心思想是：

1. **Multiplex（多重结构）**：给定一组网络（共享同一顶点集），如果取其中任意 k 个网络，可以构造：
   - **上层（up-layer）**：每条边取这 k 个网络中的最大值（加权网络）或并集（无权网络）
   - **下层（down-layer）**：每条边取这 k 个网络中的最小值或交集
   - 一个网络如果落在上、下层之间，就是"有代表性的"

2. **Edgewise 版本**：全局版本（gMD）要求网络所有边都落在 multiplex 内才计数，太严格。eMD 放宽为**逐边判断**——每条边是否落在该边的 edgewise multiplex 内，然后对所有边做加权平均，最终得到一个 0~1 的深度值。

3. **用途**：排序、找出代表性的"中心网络"、检测异常网络。

---

## 参考文献

Zhang, M., Deng, L., & Dai, W. (2025). Multiplex Depth for Network-Valued Data and Applications. *JCGS*, 34(4), 1625–1641.
