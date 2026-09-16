# -*- coding: utf-8 -*-
"""
把增量式算法嵌入连续时间 ODE 并做稳定性分析 —— 最小完整流水线
==============================================================

对任意增量式（迭代 / 在线 / 递归 / 随机逼近）算法

        x_{k+1} = x_k + eta_k * H(x_k, xi_k)          (eta_k > 0, xi_k 是数据/噪声)

统一的七步流程：

  (1) 取条件期望        h(x) = E[ H(x, xi) ]            -> 平均场方向场
  (2) 时间重标定        t = sum_i eta_i, 令 eta -> 0     -> ODE:  x' = h(x)
  (3) 平衡点            h(x*) = 0
  (4) 线性化            J = Dh(x*),  Hurwitz 判据        -> 局部指数稳定 + 收敛速率
  (5) Lyapunov 证书     A'P + P A = -Q, P > 0            -> 全局 / 非线性的稳定证明
  (6) 离散-连续桥接     x_{k+1} = (I + eta J) x_k        -> 步长上界 eta_max
                        (算法本身就是 ODE 的一个数值积分格式；见"绝对稳定域")
  (7) 噪声的影响        渐近方差 O(eta)                  -> 常数步长停在邻域, 衰减步长收敛

三个例子把这七步逐项数值验证：

  A. 线性增量迭代 (Richardson / 二次函数梯度下降)
        ODE 无步长限制，离散有 eta < 2/lambda_max
  B. Heavy-ball / momentum
        连续时间只要 beta<1 就稳定；离散还要 eta*lambda_max < 2(1+beta)
        —— 演示"连续稳定 != 离散稳定"这条最重要的边界
  C. LMS / 随机逼近
        噪声不改变 ODE 的平衡点，只决定 O(eta) 的邻域半径

运行:  py ode_incremental_stability.py
输出:  控制台表格 + figs/*.png
"""

import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
# matplotlib 字体缓存的默认位置可能不可写，指到工作区内
os.environ.setdefault("MPLCONFIGDIR", str(HERE / ".mplcache"))

import numpy as np
from scipy.linalg import expm, solve_lyapunov

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIGS = HERE / "figs"
FIGS.mkdir(exist_ok=True)
RNG = np.random.default_rng(20240607)

plt.rcParams.update({
    "figure.dpi": 130, "font.size": 9, "axes.grid": True,
    "grid.alpha": 0.3, "savefig.bbox": "tight", "legend.frameon": False,
})

RULE = "=" * 78


def spd(n, lo, hi, rng):
    """随机对称正定矩阵及其特征值（当作二次型 Hessian）。"""
    Q, _ = np.linalg.qr(rng.normal(size=(n, n)))
    lam = np.linspace(lo, hi, n)
    A = (Q * lam) @ Q.T
    return 0.5 * (A + A.T), lam


def bisect(pred, lo, hi, iters=60):
    """找 pred 由 True 变 False 的临界点。"""
    assert pred(lo) and not pred(hi)
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if pred(mid):
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ============================================================================
# A. 线性增量迭代：  x_{k+1} = x_k - eta (A x_k - b)
# ============================================================================
def demo_A():
    print(RULE)
    print("A. Linear incremental iteration  x_{k+1} = x_k - eta (A x_k - b)")
    print("   <=> gradient descent on f(x) = 1/2 x'Ax - b'x")
    print(RULE)

    n = 6
    A, _ = spd(n, 0.5, 4.0, RNG)
    lam = np.linalg.eigvalsh(A)
    b = RNG.normal(size=n)
    xstar = np.linalg.solve(A, b)
    x0 = RNG.normal(size=n)
    e0 = x0 - xstar

    # ---- (1)(2) 方向场 -> ODE ------------------------------------------------
    #   H(x) = -(Ax - b)  --期望(这里确定)-->  h(x) = -(Ax - b)
    #   ODE:  x' = -A(x - x*)       平移后: e' = -A e
    # ---- (4) 线性化 ----------------------------------------------------------
    J = -A
    hurwitz = bool(np.all(np.linalg.eigvals(J).real < 0))
    # ---- (5) Lyapunov 证书 ---------------------------------------------------
    P = solve_lyapunov(J.T, -np.eye(n))            # J'P + PJ = -I
    lyap_ok = bool(np.all(np.linalg.eigvalsh(P) > 0))

    print(f"  spectrum lambda(A)          : {np.array2string(lam, precision=3)}")
    print(f"  (4) J = -A  Hurwitz?        : {hurwitz}   (all Re(lambda) < 0)")
    print(f"  (4) ODE decay rate          : {lam.min():.4f}  (slowest mode, e^-lam_min*t)")
    print(f"  (5) Lyapunov P > 0 ?        : {lyap_ok}")
    print("  => CONTINUOUS time: stable for EVERY eta > 0 (an ODE has no step size)")

    # ---- (6) 离散桥接: M(eta) = I - eta*A, 稳定 <=> rho(M) < 1 ---------------
    rho = lambda eta: float(np.max(np.abs(1.0 - eta * lam)))
    eta_max_th = 2.0 / lam.max()
    eta_max_num = bisect(lambda e: rho(e) < 1.0, 1e-12, 4.0 / lam.max())
    print(f"  (6) eta_max   theory 2/lambda_max = {eta_max_th:.8f}")
    print(f"                numeric rho(I-etaA)=1 at {eta_max_num:.8f}  "
          f"(rel.err {abs(eta_max_num - eta_max_th) / eta_max_th:.2e})")

    # ---- 仿真验证：实测衰减率 vs 理论谱半径 --------------------------------
    def sim_ratio(eta, iters=400):
        e = e0.copy()
        prev, ratios = None, []
        for k in range(iters):
            e = (np.eye(n) - eta * A) @ e
            nrm = np.linalg.norm(e)
            if prev is not None and k > iters // 2:
                ratios.append(nrm / prev)
            prev = nrm
            if not np.isfinite(nrm) or nrm > 1e12:
                return np.inf
        return float(np.mean(ratios))

    print("\n  (6) simulated growth ratio ||e_{k+1}||/||e_k||  vs  rho(I-eta A):")
    print("         eta        measured        theory")
    for eta in (0.05, 0.15, 0.25, 0.40):
        print(f"      {eta:<9.4g}  {sim_ratio(eta):>12.8f}  {rho(eta):>12.8f}")

    # ---- 边界两侧的仿真判别 --------------------------------------------------
    print("\n  (6) straddle the boundary with the real iteration (2000 iters):")
    for tag, eta in (("0.99*eta_max", 0.99 * eta_max_th), ("1.01*eta_max", 1.01 * eta_max_th)):
        e = e0.copy()
        with np.errstate(over="ignore", invalid="ignore"):
            for _ in range(2000):
                e = (np.eye(n) - eta * A) @ e
        nrm = float(np.linalg.norm(e))
        print(f"      eta = {tag:<12s} = {eta:.6f} : ||x_2000 - x*|| = {nrm:.3e}"
              f"   ({'converges' if nrm < 1e-8 else 'DIVERGES'})")

    # ---- 图 1：时间重标定后离散轨迹收敛到 ODE 流 -----------------------------
    fig, ax = plt.subplots(figsize=(5.8, 3.7))
    tmax = 12.0
    ts = np.linspace(0.0, tmax, 300)
    ax.semilogy(ts, [np.linalg.norm(expm(-A * t) @ e0) for t in ts],
                "k-", lw=2.2, label=r"ODE  $\dot e=-Ae$  (exact)")
    for eta, style in zip((0.4, 0.2, 0.05), ("C0--", "C1-.", "C3:")):
        e = e0.copy()
        tt, err = [0.0], [np.linalg.norm(e)]
        for k in range(int(tmax / eta)):
            e = (np.eye(n) - eta * A) @ e
            tt.append((k + 1) * eta)
            err.append(np.linalg.norm(e))
        ax.semilogy(tt, err, style, lw=1.4, label=fr"discrete $\eta={eta}$")
    ax.set_xlabel(r"rescaled time   $t=k\eta$")
    ax.set_ylabel(r"$\|x_k-x^\star\|$")
    ax.set_title(r"A: discrete iterates $\to$ the ODE flow as $\eta\to0$")
    ax.legend()
    fig.savefig(FIGS / "A_ode_limit.png")
    plt.close(fig)
    return lam, eta_max_th


# ============================================================================
# B. Heavy-ball / momentum
# ============================================================================
def demo_B():
    print()
    print(RULE)
    print("B. Heavy-ball  x_{k+1} = x_k - eta*grad f(x_k) + beta (x_k - x_{k-1})")
    print(RULE)

    n = 4
    A, _ = spd(n, 0.4, 2.0, RNG)
    lam = np.linalg.eigvalsh(A)
    lmax = float(lam.max())
    x0 = RNG.normal(size=n)

    eye, zero = np.eye(n), np.zeros((n, n))

    def M(eta, beta):
        return np.block([[(1.0 + beta) * eye - eta * A, -beta * eye], [eye, zero]])

    def rho(eta, beta):
        return float(np.max(np.abs(np.linalg.eigvals(M(eta, beta)))))

    print("  (2) continuous limit (t = k*sqrt(eta), beta fixed):")
    print("        x'' + a x' + grad f = 0,     a = (1-beta)/sqrt(eta)")
    print("      energy E = 1/2||x'||^2 + f(x)  =>  dE/dt = -a||x'||^2")
    print("      => ODE stable for ALL eta>0 as long as beta < 1   (a > 0)")
    print("  (6) discrete: per mode, roots of  s^2 - (1+beta-eta*lambda)s + beta = 0")
    print("      Jury/Schur:  |beta| < 1   and   0 < eta*lambda < 2(1+beta)")
    print(f"      => eta_max(beta) = 2(1+beta)/lambda_max,   lambda_max = {lmax:.4f}")

    print("\n  empirical (numerical spectral radius) vs theoretical eta_max(beta):")
    print("       beta      theory 2(1+b)/lmax     numeric rho=1        rel.err")
    for beta in (0.0, 0.2, 0.5, 0.8, 0.9):
        th = 2.0 * (1.0 + beta) / lmax
        emp = bisect(lambda e: rho(e, beta) < 1.0, 1e-12, 6.0 / lmax, iters=70)
        print(f"      {beta:<6.2f}   {th:>16.6f}     {emp:>16.6f}      "
              f"{abs(emp - th) / th:.2e}")

    # ---- 关键点：ODE 稳定但离散发散 -----------------------------------------
    beta = 0.5
    eta_gap = 1.5 * 2.0 * (1.0 + beta) / lmax          # 明确落在稳定域之外
    a_cont = (1.0 - beta) / np.sqrt(eta_gap)
    x, xp = x0.copy(), x0.copy()
    norms, blow_at = [], None
    with np.errstate(over="ignore", invalid="ignore"):
        for k in range(400):
            xn = x - eta_gap * (A @ x) + beta * (x - xp)
            xp, x = x, xn
            nrm = float(np.linalg.norm(x))
            norms.append(nrm)
            if blow_at is None and (not np.isfinite(nrm) or nrm > 1e12):
                blow_at = k + 1
                break
    print(f"\n  GAP DEMO  beta={beta}, eta={eta_gap:.4f} > eta_max="
          f"{2*(1+beta)/lmax:.4f}")
    print(f"    continuous damping a=(1-beta)/sqrt(eta) = {a_cont:.4f} > 0 -> ODE STABLE")
    print(f"    discrete simulation: ||x_20|| = {norms[min(19, len(norms)-1)]:.4g}, "
          f"||x|| blew past 1e12 at step {blow_at} -> DIVERGES")

    # ---- 实测增长比 vs 谱半径（边界内侧） -----------------------------------
    # 用几何平均 (||s_K||/||s_k0||)^(1/(K-k0)) 估计主导模态模长，收敛到谱半径；
    # 注意用 max|entry| 而不是 norm —— 状态已经衰减到 1e-270, 平方会下溢成 0
    print("\n  simulated geometric growth rate vs rho (beta=0.5, inside):")
    for frac in (0.5, 0.8, 0.95):
        eta = frac * 2.0 * (1.0 + beta) / lmax
        z, zp = x0.copy(), x0.copy()
        k0, k1 = 600, 1800
        with np.errstate(over="ignore", invalid="ignore"):
            for k in range(k1):
                zn = z - eta * (A @ z) + beta * (z - zp)
                zp, z = z, zn
                if k == k0 - 1:
                    n0 = float(np.max(np.abs(np.concatenate([z, zp]))))
        n1 = float(np.max(np.abs(np.concatenate([z, zp]))))
        meas = float(np.exp((np.log(n1) - np.log(n0)) / (k1 - k0)))
        print(f"      eta={eta:.4f}  measured {meas:.8f}   "
              f"theory rho {rho(eta, beta):.8f}   rel.err "
              f"{abs(meas - rho(eta, beta)) / rho(eta, beta):.1e}")

    # ---- 图 2：(eta, beta) 平面上的实测稳定域 -------------------------------
    n_eta, n_beta = 46, 42
    etas = np.linspace(0.0, 4.0 / lmax, n_eta)
    betas = np.linspace(0.0, 0.98, n_beta)
    peak = np.zeros((n_beta, n_eta))
    with np.errstate(over="ignore", invalid="ignore"):
        for i, beta in enumerate(betas):
            for j, eta in enumerate(etas):
                if eta == 0.0:
                    peak[i, j] = np.log10(np.linalg.norm(x0))
                    continue
                z, zp, pk = x0.copy(), x0.copy(), np.linalg.norm(x0)
                for _ in range(3000):
                    zn = z - eta * (A @ z) + beta * (z - zp)
                    if not np.all(np.isfinite(zn)) or np.max(np.abs(zn)) > 1e8:
                        pk = 1e6
                        break
                    zp, z = z, zn
                    pk = max(pk, np.linalg.norm(z))
                peak[i, j] = min(np.log10(pk), 4.0)

    fig, ax = plt.subplots(figsize=(6.2, 4.1))
    mesh = ax.pcolormesh(etas, betas, peak, cmap="RdYlBu_r", vmin=-2, vmax=2,
                         shading="auto")
    b_line = np.linspace(0.0, 1.0, 200)
    ax.plot(2.0 * (1.0 + b_line) / lmax, b_line, "k-", lw=2.0,
            label=r"theory: $\eta\lambda_{\max}=2(1+\beta)$")
    ax.axhline(1.0, color="k", ls=":", lw=1.2)
    ax.text(0.03, 0.955, r"$\beta=1$: continuous marginal (a=0)",
            transform=ax.transAxes, fontsize=8)
    ax.annotate("GAP:\nODE stable,\ndiscrete diverges",
                xy=(eta_gap, beta), xytext=(0.66, 0.34), textcoords="axes fraction",
                fontsize=8, arrowprops=dict(arrowstyle="->", lw=1.0))
    ax.set_xlabel(r"step size  $\eta$")
    ax.set_ylabel(r"momentum  $\beta$")
    ax.set_title("B: measured stability region vs the Jury/Hurwitz prediction")
    ax.set_ylim(0, 1.0)
    ax.legend(loc="lower right")
    fig.colorbar(mesh, ax=ax, label=r"$\log_{10}\max_k\|x_k\|$")
    fig.savefig(FIGS / "B_momentum_stability.png")
    plt.close(fig)


# ============================================================================
# C. LMS / 随机逼近
# ============================================================================
def demo_C():
    print()
    print(RULE)
    print("C. LMS / stochastic approximation")
    print("   theta_{k+1} = theta_k + eta_k phi_k (y_k - phi_k' theta_k),  phi~N(0,Sigma)")
    print(RULE)

    n = 5
    Sigma, _ = spd(n, 0.6, 2.0, RNG)
    lam = np.linalg.eigvalsh(Sigma)
    L = np.linalg.cholesky(Sigma)
    theta_star = RNG.normal(size=n)
    sigma2 = 0.25
    e0 = RNG.normal(size=n)

    print("  (1) h(theta) = E[phi phi'](theta*-theta) = Sigma (theta*-theta)")
    print("  (2) ODE: theta' = Sigma (theta* - theta)")
    print(f"  (4) J = -Sigma is Hurwitz, lambda_min = {lam.min():.4f}")
    print("      (persistent excitation  <=>  Sigma > 0)")
    print(f"      => theta* exponentially stable, rate {lam.min():.4f}")

    # ---- (6) 均值 = 前向 Euler：全局误差 O(eta) ------------------------------
    print("\n  (6) mean of the iterates == forward Euler discretization of the ODE:")
    for eta in (0.05, 0.01, 0.002):
        K = int(10.0 / eta)
        Mk = np.eye(n)
        worst = 0.0
        for k in range(K + 1):
            worst = max(worst, np.linalg.norm(Mk @ e0 - expm(-Sigma * k * eta) @ e0))
            if k < K:
                Mk = Mk @ (np.eye(n) - eta * Sigma)
        print(f"      eta={eta:<7.4g} max_k ||E[theta_k] - theta_ode(t_k)|| = {worst:.3e}")

    # ---- (7) 常数步长：稳态 MSD = O(eta) ------------------------------------
    # 精确稳态协方差 V 满足（Isserlis 展开高斯四阶矩）
    #   Sigma V + V Sigma = eta [ Sigma*tr(Sigma V) + 2 Sigma V Sigma + sigma^2 Sigma ]
    # 首阶渐近 V -> (eta sigma^2 / 2) I, 即 tr V ~ eta sigma^2 n / 2
    # 这是关于 V 的线性方程 -> 直接在对称矩阵空间里组装矩阵并求解。
    # （别用 V <- ... 的定点迭代：算子 Sigma*tr(Sigma .) 是非正规的，瞬态放大
    #   会把 1e-16 的舍入误差放大到 1e140，迭代几千步后数值爆炸。）
    idx = [(i, j) for i in range(n) for j in range(i, n)]

    def _basis(i, j):
        E = np.zeros((n, n))
        E[i, j] = 1.0
        if i != j:
            E[j, i] = 1.0
        return E

    def _resid(V, eta):
        """方程左端的线性算子部分 Sigma V + V Sigma - eta[Sigma tr(Sigma V) + 2 Sigma V Sigma]"""
        SV = Sigma @ V
        return SV + SV.T - eta * (Sigma * np.trace(SV) + 2.0 * (SV @ Sigma))

    def exact_msd(eta):
        m = len(idx)
        A = np.zeros((m, m))
        for c, (i, j) in enumerate(idx):
            R = _resid(_basis(i, j), eta)
            A[:, c] = [R[a, b] for (a, b) in idx]
        b = np.array([(eta * sigma2 * Sigma)[a, b] for (a, b) in idx])
        x = np.linalg.solve(A, b)
        V = np.zeros((n, n))
        for (i, j), val in zip(idx, x):
            V[i, j] = V[j, i] = val
        return float(np.trace(V))

    R = 300
    etas = np.array([0.002, 0.004, 0.008, 0.016, 0.032, 0.064, 0.128])
    msd = np.zeros_like(etas)
    msd_th = np.array([exact_msd(e) for e in etas])
    for i, eta in enumerate(etas):
        K = max(4000, int(40.0 / (eta * lam.min())))    # burn-in 要覆盖 1/(eta*lam_min)
        burn = K // 2
        Th = np.tile(e0, (R, 1))                        # (R,n) 误差 theta - theta*
        acc = 0.0
        with np.errstate(over="ignore", invalid="ignore"):
            for k in range(K):
                phi = RNG.standard_normal((R, n)) @ L.T
                eps = RNG.normal(0.0, np.sqrt(sigma2), size=R)
                Th = Th - eta * phi * (phi * Th).sum(axis=1)[:, None] + eta * phi * eps[:, None]
                if k >= burn:
                    acc += float(np.mean(np.sum(Th ** 2, axis=1)))
        msd[i] = acc / (K - burn)
    slope_all = float(np.polyfit(np.log(etas), np.log(msd_th), 1)[0])
    slope_small = float(np.polyfit(np.log(etas[:3]), np.log(msd_th[:3]), 1)[0])
    print("\n  (7) constant step: stationary MSD is Theta(eta)")
    print("         eta      measured(MC)     exact theory   1st order eta*s^2*n/2")
    for eta, v, t in zip(etas, msd, msd_th):
        print(f"      {eta:<9.4g}  {v:>12.5f}   {t:>12.5f}   {eta * sigma2 * n / 2:>14.5f}")
    print(f"      log-log slope: exact theory {slope_all:.4f} (full range), "
          f"{slope_small:.4f} (small eta -> 1.0)")

    # ---- 衰减步长 vs 常数步长 -----------------------------------------------
    c = 0.6
    x_c, x_d = e0.copy(), e0.copy()
    err_c, err_d = [], []
    for k in range(1, 6001):
        phi = RNG.standard_normal(n) @ L.T
        y = phi @ theta_star + RNG.normal(0.0, np.sqrt(sigma2))
        x_c = x_c + 0.05 * phi * (y - phi @ (theta_star + x_c))
        x_d = x_d + (c / k) * phi * (y - phi @ (theta_star + x_d))
        err_c.append(np.linalg.norm(x_c))
        err_d.append(np.linalg.norm(x_d))

    # ---- 图 3 ----------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(7.8, 3.4))
    ax = axes[0]
    ax.loglog(etas, msd, "o", ms=4, label="Monte-Carlo")
    ax.loglog(etas, msd_th, "k-", lw=1.4, label="exact stationary theory")
    ax.loglog(etas, etas * sigma2 * n / 2, "k--", lw=1.2,
              label=r"1st order $\eta\sigma^2 n/2$")
    ax.set_xlabel(r"step size $\eta$")
    ax.set_ylabel("stationary MSD")
    ax.set_title(f"C1: noise floor is $O(\\eta)$  (slope {slope_small:.2f})")
    ax.legend()

    ax = axes[1]
    ax.semilogy(np.arange(1, 6001), err_c, lw=1.0, label=r"constant $\eta=0.05$")
    ax.semilogy(np.arange(1, 6001), err_d, lw=1.0, label=r"decaying $\eta_k=0.6/k$")
    ax.set_xlabel("iteration $k$")
    ax.set_ylabel(r"$\|\theta_k-\theta^\star\|$")
    ax.set_title("C2: constant step plateaus, decaying step converges")
    ax.legend()
    fig.savefig(FIGS / "C_lms_noise.png")
    plt.close(fig)


def main():
    demo_A()
    demo_B()
    demo_C()
    print()
    print(RULE)
    print(f"figures written to: {FIGS}")
    print(RULE)


if __name__ == "__main__":
    main()
