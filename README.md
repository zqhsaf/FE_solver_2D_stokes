# FE_solver_2D_stokes

二维稳态 **Stokes 方程**的 Taylor–Hood 混合有限元求解器（MATLAB 实现），对应讲义
*Finite Elements, Chapter 6: Finite elements for 2D steady Stokes equation*（第 1–53 页）。

- 速度 `u`：二次三角元（basis_type = 202，6 个局部基）
- 压力 `p`：一次三角元（basis_type = 201，3 个局部基）
- 粘度 ν = 1，区域 Ω = [0,1] × [−0.25, 0]

## 控制方程与弱形式

$$
-\nabla\cdot \mathbb T(u,p)=f\ \text{在 }\Omega,\qquad \nabla\cdot u=0\ \text{在 }\Omega,\qquad u=g\ \text{在 }\partial\Omega,
$$

其中 $\mathbb T(u,p)=2\nu D(u)-pI$，$D(u)=\tfrac12(\nabla u+(\nabla u)^t)$。

Galerkin 弱形式：求 $u\in U_h\times U_h,\ p\in W_h$，使

$$
a(u_h,v_h)+b(v_h,p_h)=(f,v_h),\qquad b(u_h,q_h)=0,
$$

$$
a(u,v)=\int_\Omega 2\nu D(u):D(v)\,dxdy,\qquad
b(v,p)=-\int_\Omega p\,\nabla\cdot v\,dxdy.
$$

## 解析解与右端项（讲义第 50 页）

$$
u_1=x^2y^2+e^{-y},\quad
u_2=-\tfrac23xy^3+2-\pi\sin(\pi x),\quad
p=-(2-\pi\sin(\pi x))\cos(2\pi y).
$$

$$
f_1=-2\nu x^2-2\nu y^2-\nu e^{-y}+\pi^2\cos(\pi x)\cos(2\pi y),
$$

$$
f_2=4\nu xy-\nu\pi^3\sin(\pi x)+2\pi(2-\pi\sin(\pi x))\sin(2\pi y).
$$

## 使用方法

在 MATLAB 中：

```matlab
cd('D:\matlab\FE_solver_2D_stokes');
run_error_table        % 输出第 52、53 页两张收敛性误差表，并保存 txt
plot_stokes(16)        % 绘制速度大小 / 压力的数值解与解析解对比云图（高清 png）
```

单独求解某个网格：

```matlab
[sol,P,T,Pb_u,Tb_u,Pb_p,Tb_p] = FE_solver_2D_stokes(0,1,-0.25,0,16,16);
[uLinf,uL2,uH1,pLinf,pL2,pH1] = compute_error_2D(sol,P,T, ...
    Pb_u,Tb_u,202, Pb_p,Tb_p,201);
```

## 复现结果（h = 1/8, 1/16, 1/32, 1/64）

### 表 1　速度误差（讲义第 52 页，二次三角元）

| h | ‖u−u_h‖∞ | ‖u−u_h‖₀ | \|u−u_h\|₁ |
|---|---|---|---|
| 1/8  | 1.5057e-03 | 3.5713e-04 | 2.0034e-02 |
| 1/16 | 1.9059e-04 | 4.4830e-05 | 5.0136e-03 |
| 1/32 | 2.3853e-05 | 5.6105e-06 | 1.2537e-03 |
| 1/64 | 2.9789e-06 | 7.0155e-07 | 3.1344e-04 |

收敛阶：L∞/L2 ≈ O(h³)，H1 半范 ≈ O(h²)，与讲义结论一致。

### 表 2　压力误差（讲义第 53 页，一次三角元）

| h | ‖p−p_h‖∞ | ‖p−p_h‖₀ | \|p−p_h\|₁ |
|---|---|---|---|
| 1/8  | 5.9373e-02 | 1.1890e-02 | 8.7557e-01 |
| 1/16 | 1.4973e-02 | 2.9195e-03 | 4.3617e-01 |
| 1/32 | 3.7628e-03 | 7.2893e-04 | 2.1779e-01 |
| 1/64 | 9.4337e-04 | 1.8260e-04 | 1.0885e-01 |

收敛阶：L∞/L2 ≈ O(h²)，H1 半范 ≈ O(h)，与讲义结论一致。

> 与讲义表格的绝对值对照、压力归一化方式及全部代码改动诊断，见
> `代码修改说明.md`；改造前的原始版本本地备份在 `backup_original/`（已被 .gitignore 排除，不纳入远程仓库）。
