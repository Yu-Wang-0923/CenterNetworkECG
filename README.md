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

## 输出

输入一组网络 → 得到每个网络的 EMD 深度值（0~1），基于此可以：

- **中心网络**：深度值最高的网络就是最典型的代表
- **完整排序**：所有网络按深度从高到低排序，从最典型到最异常
- **异常检测**：深度值显著低于整体的样本被标记为异常
- **边重要性**：通过结构权重，识别哪些边对中心判定贡献最大
- **两样本检验**：比较两组（如病人 vs 对照组）的深度分布是否存在显著差异

---

## 支撑材料

- **论文**：`Context/EMD/EMD.tex` / `EMD.pdf` / `EMD.md`
- **中文翻译**：`Context/EMD/EMD_zh.tex` / `EMD_zh.pdf`
- **附录**：`Context/EMD/Supplementary Material.pdf`
- **R 实现与模拟代码**：`Context/EMD/Rtool_EMD/`

---

## 参考文献

Zhang, M., Deng, L., & Dai, W. (2025). Multiplex Depth for Network-Valued Data and Applications. *JCGS*, 34(4), 1625–1641.
