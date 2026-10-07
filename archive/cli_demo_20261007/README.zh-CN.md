# 海洋热浪与珊瑚礁微塑料科研代码示例

本仓库按Lake_Microplastics_Analysis_System的科研交付形式整理：研究代码、实际数据子集、命令行demo、计算输出、期望结果、安装与复现说明。早期网页版本保存在archive/web_demo_20261006，不再作为主交付。

## 安装和运行

在Python 3.11或3.12环境中：

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install -r requirements.txt
python run_demo.py
```

该命令从29,746个真实珊瑚礁网格输入开始，重新执行等面积优先区选择、MPEI变化分析和1,000次联合模型/空间块bootstrap，再重算Table S5并绘图。不是只读取现成表格重新展示。

原分析的七个数值函数、三个SSP、四个区域、三个气候成员、5°空间块和20260804随机种子保持一致。数据保留原始行序及逐列数值；校验参考表只在计算完成后用于对照。

## 结果

输出默认保存在outputs/conservation_demo，包括新点估计、12,000条bootstrap记录、区间汇总、Table S5、网格优先级选择及PNG/SVG/PDF图件。run_checks.json明确记录哈希、数值比较、有效重复数、运行环境与耗时。完整GBR分析的有效重复数为991；不能将缺失结果补成零。

快速流程检查可用 `python run_demo.py --replicates 25 --no-plots`，但它不重现完整1,000次bootstrap的最终区间。自有数据可以通过 `--input` 指定，需满足数据字典，不能据此继续宣称与原研究冻结结果一致。

## 边界

此demo完整复现“从已归档网格派生输入到保护配置结果”的模块。原海洋模式、环境栅格处理、MaxEnt训练和完整论文仍有其他输入需求；它们的较大分析脚本保存在research_workflow，并注明外部依赖。

主线仍使用9月6日修复版HSI/MESS与2045–2055 36-tracer MPEI。各SSP使用同一未来MPEI场。MPEI是模型外部暴露，优先配置结果不是净化效果、内剂量或保护实施后的实测成效。全球MPEI变化区间包含零。

[完整说明](README.md) · [模块说明](analysis/modules/conservation_demo/README.md) · [数据字典](data/README.md) · [方法和验证](docs/methods_and_validation.md)。代码继续通过[目标仓库PR #1](https://github.com/NKUHuLab/MHW2CoralEcoSys/pull/1)提交，main合并状态以GitHub为准。
