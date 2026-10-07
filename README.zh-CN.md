# 小样本MaxEnt训练与全量栅格投影

主入口：[Physiology_MaxEnt_demo.ipynb](Physiology_MaxEnt_demo.ipynb)。已执行并保留真实模型结果、表格和图件。不再以保护配置或网页作为主Demo。

[当前全量投影Notebook下载](https://github.com/RugkeyPro/MHW2CoralEcoSys/releases/tag/v5.0.0-full-grid-maxent-demo)。完整栅格数据独立于Git发布，约596MiB，包含26个实际GeoTIFF，下载与解压均核验SHA-256。

安装Python 3.12与Java 17+（java放入PATH），创建并激活Python虚拟环境后：

```powershell
pip install -r requirements.txt
python scripts/acquire_full_grid_inputs.py
jupyter lab Physiology_MaxEnt_demo.ipynb
```

在仓库根目录打开Notebook，选择Run All Cells。逐格展示来源核验、实际输入、空间留出、原生MaxEnt训练、基线/未来投影、生理约束、图件和报告。新结果写入 `outputs/physiology_maxent_demo/`，原项目输入保持只读。

## 内容与边界

1. Acropora实际出现点1,500条、背景点10,000条，固定种子从原SWD抽取，保留父表行号。背景不是已知缺失。
2. 原七个最终变量、LQHP特征与RM=0.5，重新训练原生MaxEnt 3.4.4。按5°空间块留出，不套用旧模型AUC。
   使用原默认500次迭代上限；运行完成不等同独立确认优化器收敛或生物学有效性。训练表记录数与原生引擎去重后的样本数分别保存。
3. GFDL-ESM4基线及SSP126/245/585的完整4338×1928原生网格，每个完整环境像元调用原生MaxEnt。静态底形/粗糙度来自原完整栅格，按完全对齐的成员范围裁取；不插值、不空间下采样。缺失保留NoData。
4. 原主线约束为投影后 `HSI × Φ^0.5`；Φ为不对称温度响应乘文石饱和度sigmoid。不是生理层嵌入训练或改写MaxEnt目标函数。参数为原结构化假设，不是本Demo拟合的实验参数。
5. 比较完整无约束/有约束GeoTIFF，逐像元回读核验，并保存连续栅格地图的PNG、SVG和PDF。摘要使用同一完整环境域的网格面积加权均值，不额外限制于礁区。范围检查不等同完整MESS。

论文原生公式已核对：manuscript_第八版20260806式(1)为HSI=HSI_MaxEnt×Φ^0.5。与当前主线投影代码一致；Acropora温度/文石饱和度生理曲线24,381个值独立重算与父代码精确一致。见docs/mainline_physiology_audit.json。原论文文件只读，未修改。

仅训练一个小样本模型，不重做四类群全量优化、正式10次bootstrap或三成员集合，也不声称已证明生理约束提高预测能力。旧保护配置Notebook保留为可选下游模块。

来源、原父脚本、许可证和哈希见[模块说明](analysis/modules/physiology_maxent_demo/README.md)。[目标仓库PR #1](https://github.com/NKUHuLab/MHW2CoralEcoSys/pull/1)仍需维护者合并。
