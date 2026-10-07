# 珊瑚礁微塑料科研Notebook

主要入口是 [MHW2CoralEcoSys_demo.ipynb](MHW2CoralEcoSys_demo.ipynb)。打开Notebook即可逐步查看实际输入、计算过程、表格、日志和图件；上传的文件保留完整执行输出，GitHub上也可直接预览。

仓库继续按Lake的科研代码交付标准保留实际数据、原计算函数、期望结果与复现说明。网页和早期根目录脚本入口分别保存在archive/web_demo_20261006及archive/cli_demo_20261007。

## 安装和运行

在Python 3.11或3.12环境中：

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install -r requirements.txt
jupyter lab MHW2CoralEcoSys_demo.ipynb
```

选择该环境的Python内核，然后点击 **Run → Run All Cells**。Notebook从29,746个真实网格开始，分8个阶段完成环境参数、输入和哈希、优先区点估计、1,000次bootstrap、区间和Table S5、参考核对、绘图、保存报告。

原分析的七个数值函数、三个SSP、四个区域、三个气候成员、5°空间块和20260804随机种子保持一致。数据保留原始行序及逐列数值；校验参考表只在计算完成后用于对照。

## 结果

输出默认保存在outputs/notebook_demo，包括新点估计、12,000条bootstrap记录、区间汇总、Table S5、网格优先级选择及PNG/SVG/PDF图件。主要表格和图件也直接显示在Notebook内。run_checks.json明确记录哈希、数值比较、有效重复数、运行环境与耗时。完整GBR分析的有效重复数为991；缺失结果保留为缺失。

需要使用自己的数据或改变参数时，修改Notebook参数单元格中的INPUT_FILE、OUTPUT、N_BOOTSTRAP、RANDOM_SEED。缩短重复数或修改输入/种子会明确标记为非完整冻结设计，不继续宣称复现原最终区间。

需要自动验证Notebook时，维护者可执行 `python scripts/execute_notebook.py`；该工具逐格执行同一个.ipynb，CI也使用这一入口。

## 边界

此demo完整复现“从已归档网格派生输入到保护配置结果”的模块。原海洋模式、环境栅格处理、MaxEnt训练和完整论文仍有其他输入需求；它们的较大分析脚本保存在research_workflow，并注明外部依赖。

主线仍使用9月6日修复版HSI/MESS与2045–2055 36-tracer MPEI。各SSP使用同一未来MPEI场。MPEI是模型外部暴露，优先配置结果不是净化效果、内剂量或保护实施后的实测成效。全球MPEI变化区间包含零。

[完整说明](README.md) · [模块说明](analysis/modules/conservation_demo/README.md) · [数据字典](data/README.md) · [方法和验证](docs/methods_and_validation.md)。代码继续通过[目标仓库PR #1](https://github.com/NKUHuLab/MHW2CoralEcoSys/pull/1)提交，main合并状态以GitHub为准。
