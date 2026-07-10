# 中心网 (CenterNetworkECG)

**中心网**是一个专注于**网络型数据（network-valued data）的统计深度分析**的开源项目，旨在从一组网络中提取出最具代表性的**中心网络**，并实现网络的排序、异常检测与假设检验。

**GitHub:** [Yu-Wang-0923/CenterNetworkECG](https://github.com/Yu-Wang-0923/CenterNetworkECG)

---

## 背景：什么是 EMD？

EMD（Edgewise Multiplex Depth，边级多重深度）是 Zhang, Deng & Dai (2025) 发表于 *Journal of Computational and Graphical Statistics* 的一种针对网络型数据的**统计深度（statistical depth）**方法。它把经典的多变量 Tukey 深度推广到了网络数据领域。

### 核心思想

给定一组共享同一顶点集的网络（例如：一组相关性网络、脑功能连接网络或社交网络），EMD 通过以下方式衡量每个网络的代表性：

1. **构造 Multiplex（多重结构）**：从样本中任取 k 个网络，对每条边取它们的最大值与最小值（加权网络）/ 并集与交集（无权网络），形成上下两层。落在上下层之间的网络被认为是"有代表性的"。

2. **边级松弛（Edgewise）**：全局版本（gMD）要求网络所有边同时落在 multiplex 内才计数，对大型网络过于严格。EMD 逐边判断，将各边的结果做加权平均，得到 [0,1] 之间的深度值。

3. **结构权重**：可引入顶点中心性、边中心性等先验知识，让重要结构获得更高权重。

### 理论性质

文章证明了 EMD 满足统计深度的经典性质：

- **中心最大化**：深度值在真实的中心网络处唯一达到最大
- **单调性**：离中心越远，深度值单调递减
- **仿射不变性**
- **无穷远处消失**

并建立了样本深度估计量与深度中心的**强相合性**（Theorem 2）：当样本量增大时，EMD 值最高的网络几乎必然收敛到真实的中心网络。

### 应用场景

- 给定一组功能连接网络 → 返回最"典型"的中心网络
- 检测异常网络（深度值显著低于整体的样本）
- 两组网络之间的两样本检验
- 为网络数据提供排序与探索性分析

---

## 项目结构

```
中心网/
├── README.md
└── Context/
    └── EMD/
        ├── EMD.md               # 原始论文 Markdown
        ├── EMD.tex              # 原始论文 LaTeX
        ├── EMD.pdf              # 原始论文 PDF
        ├── EMD_zh.tex           # 中文翻译 LaTeX（可编译为 PDF）
        ├── EMD_zh.pdf           # 中文翻译 PDF
        ├── EMD.html             # 论文 HTML
        ├── images/              # 论文图表
        ├── EMD_content_list.json
        ├── EMD_content_list_v2.json
        ├── EMD_model.json
        ├── block_list.json
        └── layout.json
```

---

## 快速理解

> **给你一堆网络（比如 100 个人的相关性矩阵），EMD 能算出每个网络的"深度值"，深度值最高的那个就是最"中心"的代表网络，深度值最低的就是最"异常"的。**

---

## 参考文献

Zhang, M., Deng, L., & Dai, W. (2025). Multiplex Depth for Network-Valued Data and Applications. *Journal of Computational and Graphical Statistics*, 34(4), 1625-1641. DOI: [10.1080/10618600.2025.2475137](https://doi.org/10.1080/10618600.2025.2475137)

---

## 许可

MIT License
