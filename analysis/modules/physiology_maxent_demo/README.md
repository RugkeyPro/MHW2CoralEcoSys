# 生理约束 MaxEnt 训练与投影 Demo

主入口为仓库根目录 `Physiology_MaxEnt_demo.ipynb`，使用原生 Java MaxEnt 3.4.4 重新训练，而非载入冻结HSI结果。

实际 Acropora SWD 固定种子子集：1,500个出现记录、10,000个背景记录。按5°空间块留出约20%的块；原主线的最终七个变量、FC=LQHP、RM=0.5保留。这里训练一个小样本模型，不重新优化参数，不声称复现原正式10次bootstrap或正式模型性能。

当前采用完整4338×1928原生网格，GFDL-ESM4的baseline与SSP126/245/585每个完整环境像元都调用原生density.Project；不空间抽样、不插值点预测。静态底形/粗糙度来自原完整GeoTIFF，与成员完全对齐后裁去范围外的6行，不重采样。原2,000位置版本已移至archive/point_projection_20261007。

原主线约束是投影后 `HSI × Φ^0.5`，Φ为不对称温度高斯响应与文石饱和度sigmoid乘积；不是生理层作为训练变量，也不是修改MaxEnt目标函数。参数为原项目结构化假设，不能当作本Demo估计出的实验参数。

论文原生式(1)与此组合规则相同，见docs/mainline_physiology_audit.json。24,381个Acropora主线响应曲线值独立重算精确一致。完整源栅格独立数据包full-grid-inputs-v1包含26个GeoTIFF，每个文件均有SHA-256。源栅格不复制入Git；运行前执行scripts/acquire_full_grid_inputs.py。

需要Python requirements环境及Java 17+。`vendor/maxent.jar`来自作者维护的mrmaxent/Maxent官方仓库3.4.4归档版，许可证随包保留。原正式训练、投影R脚本按字节复制在parent_sources，仅作溯源，其原外部路径不适合直接执行。

结果、模型、日志和新图件写入 `outputs/physiology_maxent_demo/`。full_grid下每个情景保存无约束/有约束HSI、生理乘子、超训练范围变量数和原生clamping GeoTIFF。clamping是原生引擎诊断，不等同完整MESS。派生CSV中间文件只在完整数值核对并转成GeoTIFF后清理；主线原件只读。
