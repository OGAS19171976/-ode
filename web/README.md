# web/ — 交互式教学页

《连续时间稳定 ≠ 离散稳定》的四个可交互数值演示。**纯前端、零依赖、零构建**：
不装任何东西，双击 `index.html` 就能看（画布用 Canvas 2D，公式用 Unicode，不依赖 CDN，
断网也能用）。

## 目录

| 文件 | 作用 |
|---|---|
| `index.html` | 页面结构：① 时间重标定 ② 稳定域圆盘 ③ Hopf 分岔 ④ 噪声地板 ⑤ 页内自检 |
| `style.css` | 深色主题、响应式网格 |
| `stability-core.js` | **纯数值内核**（无 DOM）。浏览器挂 `globalThis.StabilityCore`，Node 里 `require` 同样可用 |
| `app.js` | 渲染与交互层：Canvas 绘图、滑块、拖拽 λ 点、读数表 |
| `test-core.js` | 对内核的 25 项断言，输出与 README §6 同源的数值证据表 |
| `smoke-dom.js` | 用最小 DOM 桩跑完整个渲染路径（17 项），验证交互不炸 |

## 跑起来

```bash
# 看页面：直接双击 index.html（或用任意静态服务器）

# 数值证据表（25 项断言）
node test-core.js

# 无浏览器冒烟测试（渲染路径 + 交互）
node smoke-dom.js
```

预期：`test-core.js` → `25 passed, 0 failed`；`smoke-dom.js` → `17 passed, 0 failed`。

> `smoke-dom.js` 只证明渲染路径不抛异常、读数写进了 DOM。
> **布局与视觉仍需在真实浏览器里确认一次**（本仓库作者未在 CI 里跑无头浏览器）。

## 部署到 GitHub Pages

页面没有任何构建步骤，把它当作静态站点根即可，三种做法任选：

1. **docs 目录法（最省事）**：把本目录内容复制到仓库的 `docs/`，然后在
   Settings → Pages → Source 选 `Deploy from a branch`，分支 `main`，目录 `/docs`。
   页面地址形如 `https://<user>.github.io/<repo>/`。
2. **分支法**：把 `web/` 内容推到 `gh-pages` 分支，Pages 源选该分支根目录。
3. **Actions 法**：加一个 `actions/upload-pages-artifact` workflow，`path: incremental-ode/web`。

三者都不需要 `npm install`、不需要改任何路径（页面内全部是相对引用）。

## 数值内核 API（`StabilityCore`）

| 分组 | 函数 | 说明 |
|---|---|---|
| 线代 | `eig2x2(m)` / `spectralRadius` / `isHurwitz` | 2×2 解析特征值与谱半径 |
| Lyapunov | `solveLyapunov2x2(J, Q)` | 解 `JᵀP + PJ = −Q`，返回 `{P, pd, residual}` |
| 前向 Euler | `eulerRho(lams, η)` / `eulerEtaMax(lams)` | 谱半径与步长上界 `min 2(−Reλ)/\|λ\|²` |
| 前向 Euler | `bisectEtaMax(rhoFn, lo, hi)` | 数值求 `ρ(η)=1` 的临界点（与解析式互验） |
| 模拟 | `simulateDiagonal` / `simulateHeavyBall` / `exactFlowDiagonal` | 离散迭代与 ODE 精确流 |
| Heavy-ball | `heavyBallEigs/ Rho / EtaMax / Disc`，`continuousDamping` | `η_max = 2(1+β)/λ_max`，复根区 `ρ=√β`，阻尼 `a=(1−β)/√η` |
| 噪声 | `steadyVariance(η, s², σ²)` = `ησ²/(2−ησ²)` | 常步长 LMS 的稳态方差（噪声地板） |
| 噪声 | `simulateLMS` / `windowStats` / `makeRng` | 可复现蒙特卡洛与标准误 |

内核被刻意写成"零 DOM"，因为它就是下一阶段 `stability-lens` 的算法原型：
同一套公式换一个封装（Python 包 + CLI），网页与库共享一份口径。

## 口径与边界

演示里的"离散 ⟹ 连续"结论都默认 Robbins–Monro / Kushner–Clark 条件成立：
`h(x)=E[H(x,ξ)]` 良定义且 Lipschitz、`Ση=∞` 且 `Ση²<∞`、ODE 有紧的渐近稳定吸引集。
以及一条最容易被忽略的：**步长上界来自积分格式的绝对稳定域，不是来自 ODE。**
