# 中心网 (CenterNetworkECG)

从一组网络（如相关性网络）中提取最具代表性的**中心网络**。

**GitHub:** [Yu-Wang-0923/CenterNetworkECG](https://github.com/Yu-Wang-0923/CenterNetworkECG)

---

## 方法

基于 EMD（Edgewise Multiplex Depth, Zhang, Deng & Dai 2025），一种网络型数据的统计深度。核心思路：

1. 对样本中任意 k 个网络，每条边取最大/最小值构成上下层（multiplex）
2. 逐边判断候选网络是否落在上下层之间，加权平均得深度值
3. **深度值最高的网络即中心网络**

详见 `Context/EMD/`。

---

## 参考文献

Zhang, M., Deng, L., & Dai, W. (2025). Multiplex Depth for Network-Valued Data and Applications. *JCGS*, 34(4), 1625–1641. DOI: [10.1080/10618600.2025.2475137](https://doi.org/10.1080/10618600.2025.2475137)
