# 生理约束 MaxEnt 训练与投影 Demo

主入口为仓库根目录 `Physiology_MaxEnt_demo.ipynb`，使用原生 Java MaxEnt 3.4.4 重新训练，而非载入冻结HSI结果。

实际 Acropora SWD 固定种子子集：1,500个出现记录、10,000个背景记录。按5°空间块留出约20%的块；原主线的最终七个变量、FC=LQHP、RM=0.5保留。这里训练一个小样本模型，不重新优化参数，不声称复现原正式10次bootstrap或正式模型性能。

2,000个固定全球背景位置，使用GFDL-ESM4的baseline、SSP126/245/585实际环境栅格抽样。保留父表行号与缺失环境值。静态bathymetry/rugosity取父SWD；其他环境量来自对应情景栅格，单位与父输入一致。

原主线约束是投影后 `HSI × Φ^0.5`，Φ为不对称温度高斯响应与文石饱和度sigmoid乘积；不是生理层作为训练变量，也不是修改MaxEnt目标函数。参数为原项目结构化假设，不能当作本Demo估计出的实验参数。

需要Python requirements环境及Java 17+。`vendor/maxent.jar`来自作者维护的mrmaxent/Maxent官方仓库3.4.4归档版，许可证随包保留。原正式训练、投影R脚本按字节复制在parent_sources，仅作溯源，其原外部路径不适合直接执行。

结果、模型、日志和新图件写入 `outputs/physiology_maxent_demo/`。主线原件只读。
