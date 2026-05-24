# 仿真部分审计报告

本报告检查 reproduction 输出中的 0D 求解结果和 3D svMultiPhysics 输入文件。

结论：0D 已经实际运行并生成 `output.csv`；3D 目前已生成网格与 XML 输入文件，但本机未发现 `svMultiPhysics` 可执行程序，也未发现 3D `result*` 输出，因此当前不能声明 3D CFD 已完成。

## 求解器状态

- `svMultiPhysics`: `not found in PATH`
- 0D solver candidates:
  - `/root/autodl-tmp/svVascularize/svv/bin/svzerodsolver` exists=False
  - `/root/autodl-tmp/svVascularize/svv/utils/solvers/0D/Linux/x86_64/svzerodsolver` exists=True
  - `/root/autodl-tmp/svVascularize/svv/solvers/svzerodsolver` exists=False

## 0D 求解结果

| 体系 | rows | vessels | timepoints | flow balance rel. error | pressure_in range | pressure_out range |
|---|---:|---:|---:|---:|---:|---:|
| `cube` | 35 | 7 | 5 | 4.947e-10 | 21938.2 to 50310.3 | -5.0822e-21 to 42544.6 |
| `cube_3d` | 15 | 3 | 5 | 3.590e-10 | 38332 to 50310.3 | -5.0822e-21 to 38332 |
| `cube_forest` | 220 | 44 | 5 | 8.056e-10 | 25323.2 to 147557 | -1.27055e-21 to 141185 |
| `biventricle` | 15 | 3 | 5 | 4.557e-08 | 48359 to 50310.3 | -3.04932e-20 to 48359 |

## 3D 输入文件

| 体系 | mesh exists | steps | dt | faces | missing BC faces | result files |
|---|---:|---:|---:|---:|---|---:|
| `cube` | True | 1 | 0.001 | 6 | none | 0 |
| `cube_3d` | True | 1 | 0.001 | 4 | none | 0 |
| `cube_forest` | True | 1 | 0.001 | 3 | none | 0 |
| `biventricle` | True | 1 | 0.001 | 4 | none | 0 |

## 判断

- 0D：已经完成一次 steady/准稳态 lumped-parameter 求解；这些结果可以用于检查网络阻力、压力和流量分配是否数量级合理。
- 3D：当前处于 simulation input/export 阶段；XML、体网格和面文件已经存在，但没有本机求解器和结果文件。
- 当前 XML 已补齐 `lumen_*` 血管壁面的 no-slip Dirichlet 边界条件；后续重新导出也会由代码自动补齐。
- 下一步若要做真正 3D CFD，需要安装或提供 `svMultiPhysics`，然后从最小的 `cube_3d` 开始试跑 1-10 step，检查收敛、质量守恒、压力/速度场，再扩大到 `cube_forest` 和 `biventricle`。
