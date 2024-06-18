# 吉他物理仿真

## 记号

| 记号          | 意义                                    | 值                                                           |
| ------------- | --------------------------------------- | ------------------------------------------------------------ |
| $ u_s(x, t) $ | $ t $ 时刻，弦在位置 $ x $ 处的法向位移 |                                                              |
| $ \rho_s $    | 弦的线密度                              | $0.00525$ $ \text{kg}\cdot\text{m}^{-1}$                     |
| $T$           | 弦的张力                                | $60$ $\text{N}$                                              |
| $ \eta_s $    | 黏弹性阻尼系数                          | $9\cdot10^{-8}$ $\text{s}$                                   |
| $ R_s $       | 流体阻尼系数                            | $0.75$ $\text{s}^{-1}$                                       |
| $ f_s(x, t) $ | 激励力                                  | $g(x)h(t)$                                                   |
| $g(x)$        | 激励力（位置分量）                      | $\displaystyle\frac{e^{-(\frac{x-x_0}{\delta_s})^2}}{\int_0^{l_s}e^{-(\frac{x-x_0}{\delta_s})^2}}$ |
| $h(t)$        | 激励力（时间分量）                      | $\displaystyle \begin{cases}1-\cos(\pi t/t_1), & 0\le t\le 1 \\1+\cos(\pi(t-t_1)/(t_2-t_1)), & t_1\le t\le t_2 \\ 0, & t>t_2 \end{cases}$ |
| $x_0$         | 拨弦位置                                | $55$ $\text{cm}$                                             |
| $\delta_c$    | 激励力扩散常数                          | $0.006$ $\text{m}$                                           |
| $t_1$         | 激励最强时刻                            | $0.0004$ $\text{s}$                                          |
| $t_2$         | 激励终止时刻                            | $0.015$ $\text{s}$                                           |


| 域 | 物理意义 | 定义 | 空间 | 离散化 |
| -- | ------- | --- | --- | ----- |
| $l$ | 弦 | $[0, l_s]$ | $\mathbb{R}$ | 线段，顶点 |
| $\omega$ | 前板 | | $\mathbb{R}^2$ | 三角形网格，顶点、中心、边中点 |
| $\gamma_0$ | 前板外轮廓 | | $\mathbb{R}$ | |
| $\gamma_f$ | 前板内轮廓 | | $\mathbb{R}$ | |
| $\Sigma$ | 侧板 | | $\mathbb{R}^2$ | 
| $\Gamma$ | 板 | $\omega \cup \Sigma$ | $\mathbb{R}^2$ |
| $\Omega$ | 空气 | $\mathbb{R}^3 \setminus \Gamma$ | $\mathbb{R}^3$ | 立方体网格，中心 |


| 未知量 | 物理意义 | 域映射 | 定义 | 有限元 | 系数向量 |
| ----- | ------- | ----- | --- | ----- | ------- |
| $v_s$   | 弦的法向速度 | $\mathbb{R}\rightarrow\mathbb{R}$ | $\partial_t u_s$ | $P_0$ (piecewise constant) | $v_{s_h}$ |
| $q$    | 弦的张力 | $\mathbb{R}\rightarrow\mathbb{R}$ | $T\partial_x u_s$ | $P_1$ (piecewise linear) |
| $v_p$ | 板的法向速度 | $\omega\rightarrow\mathbb{R}$ | $\partial_t u_p$ | $P_2$ (lagrange) |
| $M$ | 板的弯曲矩 | $\omega\rightarrow\mathbb{R}^3$ | $a^3 \mathbf{C}\varepsilon(\nabla u_p)$ | $P_2$ (lagrange) |
| $\lambda$ | 接触面压力差 | $\Gamma\rightarrow\mathbb{R}$ | $p_e - p_i$ | $P_0$ (lagrange) |
| $p$ | 声场的声压 | $\Omega\rightarrow\mathbb{R}$ | | $P_0$ (piecewise constant in cube) |
| $\mathbf{v}_a$ | 声场的空气速度 | $\Omega\rightarrow\mathbb{R}^3$ | | $P_0$ (Raviart–Thomas in cube)|

| 有限元系数向量 | 物理意义 | 域映射 | 定义 | 有限元 |
| ----- | ------- | ----- | --- | ----- |
| $v_{s_h}$   | 弦的法向速度 | $\mathbb{R}\rightarrow\mathbb{R}$ | $\partial_t u_s$ | $P_0$ (piecewise constant) |
| $q_h$    | 弦的张力 | $\mathbb{R}\rightarrow\mathbb{R}$ | $T\partial_x u_s$ | $P_1$ (piecewise linear) |
| $v_{p_h}$ | 板的法向速度 | $\omega\rightarrow\mathbb{R}$ | $\partial_t u_p$ | $P_2$ (lagrange) |
| $M_h$ | 板的弯曲矩 | $\omega\rightarrow\mathbb{R}^3$ | $a^3 \mathbf{C}\varepsilon(\nabla u_p)$ | $P_2$ (lagrange) |
| $\lambda_h$ | 接触面压力差 | $\Gamma\rightarrow\mathbb{R}$ | $p_e - p_i$ | $P_0$ (lagrange) |
| $p_h$ | 声场的声压 | $\Omega\rightarrow\mathbb{R}$ | | $P_0$ (piecewise constant in cube) |
| $\mathbf{v}_{a_h}$ | 声场的空气速度 | $\Omega\rightarrow\mathbb{R}^3$ | | $P_0$ (Raviart–Thomas in cube)|

## 模型

### 1. 弦模型
  $$
    \begin{aligned}
      & 
      \rho_s \frac{\partial^2 u_s}{\partial t^2} - T \left( 1 + \eta_s \frac{\partial}{\partial t} \right) \frac{\partial^2 u_s}{\partial x^2} + \rho_s R_s \frac{\partial u_s}{\partial t} = f_s(x, t) 
      \\ \Rightarrow \quad & 
      \begin{cases}
        \rho_s \partial_t v_s - \partial_x q = f_s \\
        \partial_t q - T \partial_x v_s = 0
      \end{cases}
      \\ \Rightarrow \quad &
      \begin{cases}
        \frac{d}{dt}\int_0^{l_s}\rho_s v_s v_s^* \, dx - 
        \int_0^{l_s} \partial_x q v_s^* = \int_0^{l_s} f_s v_s^* \, dx \\
        \frac{d}{dt}\int_0^{l_s}\frac{1}{T}qq^* + \int_0^{l_s} \partial_x q^* v_s - q^*(l_s, t) v_p(x_0, y_0) = 0
      \end{cases}
      \\ \Rightarrow \quad &
      \begin{cases}
        M_h^q \frac{d q_h}{dt} - D_h q_h = f_{s_h} \\
        M_h^v \frac{d v_{sh}}{dt} - D_h^T v_{sh} + J_h^T v_{ph} = 0
      \end{cases} , \quad \text{where} \quad \begin{cases}
        M_h^q[i][j] = \int_0^{l_s} N_i N_j \, dx \\
        D_h[i][j] = 
      \end{cases}
    \end{aligned}
  $$

### 2. 共鸣板模型

  $$
  a \rho_p \frac{\partial^2 u_p}{\partial t^2} + \left( 1 + h_p \frac{\partial}{\partial t} \right)  \nabla \cdot \nabla a^3 \mathbf{C} \varepsilon(\nabla u_p) + a r_p R_p \frac{\partial u_p}{\partial t} = F - [p]_\omega, \quad \text{in} \quad \omega
  $$
  其中：
  - $ u_p(x, y, t) $：共鸣板在位置 $ (x, y) $ 处、时间 $ t $ 的横向位移
  - $ a $：共鸣板厚度
  - $ r_p $：共鸣板密度
  - $ h_p $、$ R_p $：阻尼系数
  - $ F $：弦对共鸣板的作用力
  - $ [p]_\nu $：声场对共鸣板的压力差

### 3. 声场模型
#### 微分方程
- **方程**：
  $$
  \frac{\partial p}{\partial t} = -c_a^2 \rho_a \nabla \cdot \mathbf{v}_a
  $$
  $$
  \rho_a \frac{\partial \mathbf{v}_a}{\partial t} = -\nabla p
  $$
  其中：
  - $ p(x, y, z, t) $：声压
  - $ \mathbf{v}_a(x, y, z, t) $：声速
  - $ c_a $：空气中的声速
  - $ \rho_a $：空气密度

- **边界条件**：
  $$
  \mathbf{v}_a \cdot \mathbf{e}_z = \frac{\partial u_p}{\partial t} \quad \text{在} \quad \nu
  $$

- **初始条件**：
  $$
  p(x, y, z, 0) = p^0(x, y, z), \quad \mathbf{v}_a(x, y, z, 0) = \mathbf{v}_a^0(x, y, z)
  $$

#### 离散化版本
- **空间离散化**：
  - 使用有限差分法进行空间离散化，采用分段常数（P0）基函数 $ \chi_k(x, y, z) $ 表示压力 $ p $。

  $$
  p^h(x, y, z, t) = \sum_{k=1}^{P} p^k(t) \chi_k(x, y, z)
  $$

- **时间离散化**：
  - 使用显式中心差分法进行时间离散化：

  $$
  \frac{p^{n+1} - 2p^n + p^{n-1}}{\Delta t^2} = -c_a^2 \rho_a \nabla \cdot \mathbf{v}_a^{n}
  $$
  
  $$
  \rho_a \frac{\mathbf{v}_a^{n+1} - \mathbf{v}_a^{n-1}}{2\Delta t} = -\nabla p^{n}
  $$

## 有限元


原未知量 $(u_s, u_p, p, \mathbf{v}_a)$

| 未知量 | 物理意义 | 域映射 |
| ----- | ------- | ----- |
| $u_s$   | 弦的法向位移 | $l\rightarrow\mathbb{R}$ |
| $u_p$ | 板的法向位移 | $\omega\rightarrow\mathbb{R}$ |
| $p$ | 声场的声压 | $\Omega\rightarrow\mathbb{R}$ |
| $\mathbf{v}_a$ | 声场的空气速度 | $\Omega\rightarrow\mathbb{R}^3$ | 


等效未知量 $(v_s, q, v_p, M, \lambda, p, \mathbf{v}_a)$


## 变分形式和矩阵形式对照表

### 共鸣板模型
$$\int_{\Omega} a \rho_p v_p \frac{\partial v_p^*}{\partial t} \, d\Omega - \int_{\Omega} (\text{Div} M) \cdot (\nabla v_p^*) \, d\Omega + \int_{\Omega} a \rho_p R_p v_p v_p^* \, d\Omega = \int_{\Omega} (F - [p]_\omega) v_p^* \, d\Omega$$
$$M_p \frac{d v_{ph}}{dt} - H^T M_h + R_p M_p v_{ph} = -J_h q_h - (B_v)^T \lambda_h$$
$$\begin{aligned}
  M_p&: N_p \times N_p \\
  v_{ph}&: N_p \times 1 \\
  H^T&: N_p \times (3N_m) \\
  M_h&: 3N_m \times 1 \\
  R_p&: N_p \times N_p \\
  J_h&: N_p \times N_s \\
  q_h&: N_s \times 1 \\
  B_v&: N_p \times N_{\lambda} \\
  \lambda_h&: N_{\lambda} \times 1
\end{aligned}$$

其中
$$
M = \begin{bmatrix}
M_{xx} & M_{xy} \\
M_{xy} & M_{yy}
\end{bmatrix}
$$
$$
\text{Div} M = \left( \frac{\partial M_{xx}}{\partial x} + \frac{\partial M_{xy}}{\partial y}, \frac{\partial M_{xy}}{\partial x} + \frac{\partial M_{yy}}{\partial y} \right)^T
$$
$$
\int_{\Omega} (\text{Div} M) \cdot (\nabla v_p^*) \, d\Omega = \int_{\Omega} \left( \frac{\partial M_{xx}}{\partial x} + \frac{\partial M_{xy}}{\partial y} \right) \frac{\partial v_p^*}{\partial x} + \left( \frac{\partial M_{xy}}{\partial x} + \frac{\partial M_{yy}}{\partial y} \right) \frac{\partial v_p^*}{\partial y} \, d\Omega
$$

矩阵形式

$$
\mathbf{K} \mathbf{v_p^*} = \mathbf{f}
$$

其中 $\mathbf{K}$ 是刚度矩阵，$\mathbf{v_p^*}$ 是节点值的向量，$\mathbf{f}$ 是载荷向量。对于 $\text{Div} M$ 和 $\nabla v_p^*$ 的点乘，刚度矩阵 $\mathbf{K}$ 的元素可以表示为：

$$
K_{ij} = \int_{\Omega} \left( \frac{\partial M_{xx}}{\partial x} + \frac{\partial M_{xy}}{\partial y} \right) \frac{\partial N_i}{\partial x} + \left( \frac{\partial M_{xy}}{\partial x} + \frac{\partial M_{yy}}{\partial y} \right) \frac{\partial N_i}{\partial y} \, d\Omega
$$




### 共鸣板模型
$$\int_{\Omega} a^{-3} C_e^{-1} M \frac{\partial M^*}{\partial t} \, d\Omega + \int_{\Omega} (\text{Div} M^*) \cdot (\nabla v_p) \, d\Omega = 0$$
$$M_M \frac{d M_h}{dt} + H v_{ph} = 0$$
$M_M: 3N_m \times 3N_m$, $\frac{d M_h}{dt}: 3N_m \times 1$, $H: 3N_m \times N_p$, $v_{ph}: N_p \times 1$

### 弦模型
$$\int_0^{l_s} \rho_s v_s \frac{\partial v_s^*}{\partial t} \, dx - \int_0^{l_s} q \frac{\partial v_s^*}{\partial x} \, dx = \int_0^{l_s} f_s v_s^* \, dx$$
$M_s \frac{d v_{sh}}{dt} - D_h q_h = f_{sh}$
$M_s: N_s \times N_s$, $\frac{d v_{sh}}{dt}: N_s \times 1$, $D_h: N_s \times N_s$, $q_h: N_s \times 1$, $f_{sh}: N_s \times 1$

### 弦-板耦合
$$\int_0^{l_s} \frac{1}{T} q \frac{\partial q^*}{\partial t} \, dx + \int_0^{l_s} v_s \frac{\partial q^*}{\partial x} \, dx - q(l_s, t) v_p(x_0, y_0) = 0$$
$$M_q \frac{d q_h}{dt} + D^T v_{sh} - J^T v_{ph} = 0$$
$M_q: N_s \times N_s$, $\frac{d q_h}{dt}: N_s \times 1$, $D^T: N_s \times N_s$, $v_{sh}: N_s \times 1$, $J^T: N_s \times N_p$, $v_{ph}: N_p \times 1$

### 声场模型
$$\int_V \frac{1}{\rho_a c_a^2} p \frac{\partial p^*}{\partial t} \, dV + \int_V p^* (\text{div} v_a) \, dV = 0$$
$$M_{pa} \frac{d p_h}{dt} + G^T v_{ah} = 0$$
$M_{pa}: N_p \times N_p$, $\frac{d p_h}{dt}: N_p \times 1$, $G^T: N_p \times N_a$, $v_{ah}: N_a \times 1$

### 声场模型
$$\int_V \rho_a v_a \frac{\partial v_a^*}{\partial t} \, dV - \int_V p (\nabla \cdot v_a^*) \, dV = \int_{\Gamma} \lambda v_a^* \cdot \mathbf{n} \, d\Gamma$$
$$M_a \frac{d v_{ah}}{dt} - G_h p_h - (B_G)^T \lambda_h = 0$$
$M_a: N_a \times N_a$, $\frac{d v_{ah}}{dt}: N_a \times 1$, $G_h: N_a \times N_p$, $p_h: N_p \times 1$, $B_G^T: N_a \times N_{\lambda}$, $\lambda_h: N_{\lambda} \times 1$
$B_v v_{ph} - B_G v_{ah} = 0$
$B_v v_{ph} - B_G v_{ah} = 0$
$B_v: N_{\lambda} \times N_p$, $v_{ph}: N_p \times 1$, $B_G: N_{\lambda} \times N_a$, $v_{ah}: N_a \times 1$


