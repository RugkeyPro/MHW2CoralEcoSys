# 生理约束 MaxEnt：真实数据训练与投影 Demo

主入口：[Physiology_MaxEnt_demo.ipynb](Physiology_MaxEnt_demo.ipynb)。已执行并保留真实模型结果、表格和图件。不再以保护配置或网页作为主Demo。

[完整Notebook、实际数据与MaxEnt引擎下载](https://github.com/RugkeyPro/MHW2CoralEcoSys/releases/tag/v4.0.0-physiology-maxent-demo)。Windows与Linux完整执行均通过，13项测试通过；独立来源核验见docs/maxent_input_validation.json。

安装Python 3.12与Java 17+（java放入PATH），创建并激活Python虚拟环境后：

```powershell
pip install -r requirements.txt
jupyter lab Physiology_MaxEnt_demo.ipynb
```

在仓库根目录打开Notebook，选择Run All Cells。逐格展示来源核验、实际输入、空间留出、原生MaxEnt训练、基线/未来投影、生理约束、图件和报告。新结果写入 `outputs/physiology_maxent_demo/`，原项目输入保持只读。

## 内容与边界

1. Acropora实际出现点1,500条、背景点10,000条，固定种子从原SWD抽取，保留父表行号。背景不是已知缺失。
2. 原七个最终变量、LQHP特征与RM=0.5，重新训练原生MaxEnt 3.4.4。按5°空间块留出，不套用旧模型AUC。
3. 同一2,000个全球背景位置，抽取GFDL-ESM4基线及SSP126/245/585真实栅格环境。底形/粗糙度来自原SWD，其他量来自对应情景栅格。缺失保留为空。
4. 原主线约束为投影后 `HSI × Φ^0.5`；Φ为不对称温度响应乘文石饱和度sigmoid。不是生理层嵌入训练或改写MaxEnt目标函数。参数为原结构化假设，不是本Demo拟合的实验参数。
5. 比较无约束/有约束投影，保存PNG、可编辑SVG和PDF。抽样点不是连续全球栅格，简单均值不是全球礁区面积加权均值。范围检查不等同完整MESS。

仅训练一个小样本模型，不重做四类群全量优化、正式10次bootstrap或三成员集合，也不声称已证明生理约束提高预测能力。旧保护配置Notebook保留为可选下游模块。

来源、原父脚本、许可证和哈希见[模块说明](analysis/modules/physiology_maxent_demo/README.md)。[目标仓库PR #1](https://github.com/NKUHuLab/MHW2CoralEcoSys/pull/1)仍需维护者合并。
