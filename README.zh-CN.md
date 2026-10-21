# 二维相分离

用 Cahn–Hilliard 方程模拟二元混合物的自旋分解：从带有小扰动的初始组分场出发，观察富集相形成与畴结构粗化。

## 运行

需要 Python 3.11 或以上版本。在项目目录运行：

```powershell
py -m pip install -e .
phase-separation --output run --steps 2000
```

输出包括组分演化图、自由能及平均组分漂移图、NPZ 场数据、CSV 诊断数据和运行参数。输出目录必须为空，防止覆盖已有结果。

改变平均组分可观察非对称相分离：

```powershell
phase-separation --mean 0.3 --seed 12 --steps 4000 --output off-critical
```

## 模型与数值方法

新增案例可一次运行：

```powershell
py run_cases.py --output case-runs
```

包括对称组分、非对称组分、稳定混合物和固定终止时间的时间步对比。参数、结果与适用范围见 [案例说明](examples/CASES.md)。

## 控制方程

模型采用双阱自由能 `(phi²-1)²/4`，化学势为 `mu=phi³-phi-kappa∇²phi`，演化方程为 `∂phi/∂t=M∇²mu`。所有参数均为无量纲量，采用二维周期边界和常迁移率。

空间采用 Fourier 离散，时间采用带线性稳定项的半隐式格式。零波数模态保持不变，因此平均组分在浮点误差范围内守恒。非线性项直接在网格上计算，未做去混叠处理。程序不裁剪组分值，也不人为修改自由能曲线。

这是一份扩散驱动的相场代码，不包含 LBM 流场、温度演化或成核模型。默认参数用于演示相分离，不对应某种材料的定量预测。较大的时间步或极端初值可能失稳；定量研究需在固定域尺寸下加密网格、减小时间步并比较统计量。

## 检查结果

- `diagnostics.csv`：检查平均组分、总质量、自由能、组分极值与方差。
- `fields.npz`：`phi` 的维度为 `[保存时刻, y, x]`，`time` 为对应时间。
- `report.json`：记录保存时刻的平均组分漂移和自由能变化；保存间隔内的变化不会被检测到。
- 默认只保存每 200 步；短时精细检查可使用 `--every 1`。

```powershell
py -m unittest discover -s tests -v
```

参数表、离散更新式和测试范围见 [英文说明](README.md)。许可证为 MIT。
