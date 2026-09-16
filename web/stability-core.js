/*!
 * stability-core.js — 增量式算法 → ODE → 稳定性分析：纯数值内核
 * ---------------------------------------------------------------------------
 * 无 DOM、无依赖。浏览器 <script> 直接引入，Node 里 `require` 也一样能拿到：
 * 一切挂在 globalThis.StabilityCore 上。
 *
 * 覆盖 README 的七步流程里全部「算得出来」的部分：
 *   (3) 平衡点由调用方给出（线性化后只关心偏差 e = x - x*，故写 ẋ = J x）
 *   (4) Hurwitz 判据 / 谱半径
 *   (5) Lyapunov 方程 J'P + PJ = -Q
 *   (6) 前向 Euler 绝对稳定域 → eta_max，Heavy-ball 的 2(1+beta)/lambda_max
 *   (7) 常步长 LMS 的稳态方差（噪声地板 ∝ eta）
 */
(function (root) {
  'use strict';

  /* ======================= 复数小工具 ======================= */
  const cabs = (z) => Math.hypot(z.re, z.im);
  const cscale = (z, a) => ({ re: z.re * a, im: z.im * a });
  const cstr = (z, d) => {
    d = d === undefined ? 3 : d;
    return (z.im >= 0 ? z.re.toFixed(d) + '+' + z.im.toFixed(d) : z.re.toFixed(d) + '-' + (-z.im).toFixed(d)) + 'i';
  };

  /* ======================= 2x2 实矩阵解析特征值 ======================= */
  /** m = [[a,b],[c,d]] → [{re,im}, {re,im}] */
  function eig2x2(m) {
    const a = m[0][0], b = m[0][1], c = m[1][0], d = m[1][1];
    const tr = a + d, det = a * d - b * c;
    const disc = tr * tr - 4 * det;
    if (disc >= 0) {
      const s = Math.sqrt(disc);
      return [{ re: (tr + s) / 2, im: 0 }, { re: (tr - s) / 2, im: 0 }];
    }
    const s = Math.sqrt(-disc);
    return [{ re: tr / 2, im: s / 2 }, { re: tr / 2, im: -s / 2 }];
  }

  const spectralRadius = (eigs) => Math.max.apply(null, eigs.map(cabs));
  const isHurwitz = (eigs) => eigs.every((z) => z.re < 0);

  /* ======================= (5) Lyapunov 证书 ======================= */
  /** 解 J'P + PJ = -Q（P 对称，2x2）。返回 {P:[[p,q],[q,r]], pd:bool, residual:number} */
  function solveLyapunov2x2(J, Q) {
    const a = J[0][0], b = J[0][1], c = J[1][0], d = J[1][1];
    // 未知量 [p, q, r]
    const M = [
      [2 * a, 2 * c, 0],
      [b, a + d, c],
      [0, 2 * b, 2 * d],
    ];
    const rhs = [-Q[0][0], -Q[0][1], -Q[1][1]];
    const sol = solve3(M, rhs);
    if (!sol) return null;
    const P = [[sol[0], sol[1]], [sol[1], sol[2]]];
    const pd = sol[0] > 0 && sol[0] * sol[2] - sol[1] * sol[1] > 0;
    // 残差 ||J'P + PJ + Q||_F
    const JP = matmul(transpose(J), P), PJ = matmul(P, J);
    let res = 0;
    for (let i = 0; i < 2; i++)
      for (let j = 0; j < 2; j++) {
        const v = JP[i][j] + PJ[i][j] + Q[i][j];
        res += v * v;
      }
    return { P, pd, residual: Math.sqrt(res) };
  }

  function transpose(m) { return [[m[0][0], m[1][0]], [m[0][1], m[1][1]]]; }
  function matmul(A, B) {
    return [
      [A[0][0] * B[0][0] + A[0][1] * B[1][0], A[0][0] * B[0][1] + A[0][1] * B[1][1]],
      [A[1][0] * B[0][0] + A[1][1] * B[1][0], A[1][0] * B[0][1] + A[1][1] * B[1][1]],
    ];
  }
  /** 3x3 高斯消元（部分选主元），奇异返回 null */
  function solve3(A, y) {
    const M = [[A[0][0], A[0][1], A[0][2], y[0]], [A[1][0], A[1][1], A[1][2], y[1]], [A[2][0], A[2][1], A[2][2], y[2]]];
    for (let col = 0; col < 3; col++) {
      let piv = col;
      for (let r = col + 1; r < 3; r++) if (Math.abs(M[r][col]) > Math.abs(M[piv][col])) piv = r;
      if (Math.abs(M[piv][col]) < 1e-300) return null;
      const t = M[col]; M[col] = M[piv]; M[piv] = t;
      for (let r = 0; r < 3; r++) {
        if (r === col) continue;
        const f = M[r][col] / M[col][col];
        if (!f) continue;
        for (let k = col; k < 4; k++) M[r][k] -= f * M[col][k];
      }
    }
    return [M[0][3] / M[0][0], M[1][3] / M[1][1], M[2][3] / M[2][2]];
  }

  /* ======================= (6) 前向 Euler 绝对稳定域 ======================= */
  /** x_{k+1} = (I + eta*J) x_k，lam 为 J 的一个特征值 */
  function eulerMultiplier(lam, eta) { return { re: 1 + eta * lam.re, im: eta * lam.im }; }
  function eulerStable(lam, eta) { return cabs(eulerMultiplier(lam, eta)) < 1; }
  function eulerRho(lams, eta) { return Math.max.apply(null, lams.map((l) => cabs(eulerMultiplier(l, eta)))); }

  /**
   * 步长上界：|1+eta*lam| < 1  <=>  0 < eta < 2*(-Re lam)/|lam|^2
   * 要求 J 是 Hurwitz 的（否则连续时间本身就不稳定，返回 etaMax = 0）。
   */
  function eulerEtaMax(lams) {
    let best = Infinity, hurwitz = true;
    for (const l of lams) {
      if (!(l.re < 0)) { hurwitz = false; continue; }
      const v = 2 * (-l.re) / (l.re * l.re + l.im * l.im);
      if (v < best) best = v;
    }
    if (!hurwitz) return { etaMax: 0, hurwitz: false };
    return { etaMax: best, hurwitz: true };
  }

  /** 用二分法数值求 rho(eta)=1 的临界点（与解析 eta_max 对照用） */
  function bisectEtaMax(rhoFn, lo, hi, iters) {
    iters = iters || 60;
    if (rhoFn(lo) >= 1 || rhoFn(hi) < 1) return NaN;
    for (let i = 0; i < iters; i++) {
      const mid = 0.5 * (lo + hi);
      if (rhoFn(mid) < 1) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
  }

  /** 把离散步的有效连续衰减率还原出来：-log|1+eta*lam| / eta（小 eta 时 → -Re lam） */
  function eulerEffectiveRate(lam, eta) {
    const m = cabs(eulerMultiplier(lam, eta));
    return m > 0 ? -Math.log(m) / eta : Infinity;
  }

  /* ======================= 模拟：对角线性系统 ======================= */
  /** lams 为正的曲率（x ← (1-eta*lam)x），返回 [x_0 .. x_K] */
  function simulateDiagonal(lams, eta, steps, x0) {
    const out = [x0.slice()];
    let x = x0.slice();
    for (let k = 0; k < steps; k++) {
      for (let i = 0; i < x.length; i++) x[i] *= 1 - eta * lams[i];
      out.push(x.slice());
      for (const v of x) if (!isFinite(v) || Math.abs(v) > 1e12) return out;
    }
    return out;
  }
  /** 精确流 x_i(t) = x0_i * e^{-lam_i t} */
  function exactFlowDiagonal(lams, t, x0) { return x0.map((v, i) => v * Math.exp(-lams[i] * t)); }
  const norm = (x) => Math.sqrt(x.reduce((s, v) => s + v * v, 0));

  /* ======================= Heavy-ball / momentum ======================= */
  /** 单模态：M = [[1+beta-eta*lam, -beta], [1, 0]] */
  function heavyBallMatrix(eta, beta, lam) { return [[1 + beta - eta * lam, -beta], [1, 0]]; }
  function heavyBallEigs(eta, beta, lam) { return eig2x2(heavyBallMatrix(eta, beta, lam)); }
  function heavyBallRho(eta, beta, lam) { return spectralRadius(heavyBallEigs(eta, beta, lam)); }
  /** Jury/Schur：|beta|<1 且 0 < eta*lam < 2(1+beta) */
  function heavyBallEtaMax(beta, lamMax) { return 2 * (1 + beta) / lamMax; }
  /** 判别式 tr^2 - 4*det，<0 时两根为复共轭、|mu| = sqrt(beta) */
  function heavyBallDisc(eta, beta, lam) { const tr = 1 + beta - eta * lam; return tr * tr - 4 * beta; }
  /** 连续极限 ẍ + a ẋ + ∇f = 0 的阻尼系数；beta<1 ⟹ a>0 ⟹ ODE 无条件稳定 */
  function continuousDamping(eta, beta) { return (1 - beta) / Math.sqrt(eta); }
  /** x_{k+1} = x_k - eta*lam*x_k + beta*(x_k - x_{k-1})，x_{-1}=x_0=x0 */
  function simulateHeavyBall(eta, beta, lam, steps, x0) {
    const out = [x0];
    let xp = x0, x = x0;
    for (let k = 0; k < steps; k++) {
      const xn = x - eta * lam * x + beta * (x - xp);
      xp = x; x = xn;
      out.push(x);
      if (!isFinite(x) || Math.abs(x) > 1e12) break;
    }
    return out;
  }

  /* ======================= (7) 噪声地板：标量 LMS/SA ======================= */
  /**
   * theta_{k+1} = theta_k + eta*phi*(y - phi*theta_k),  phi^2 = s2（确定性回归子）,
   * y = phi*theta* + eps, Var[eps] = sigma2.
   * 误差递推 e_{k+1} = (1-eta*s2)e_k + eta*sqrt(s2)*eps
   * 稳态方差 V = eta^2*s2*sigma2 / (1-(1-eta*s2)^2) = eta*sigma2/(2-eta*s2)
   */
  function steadyVariance(eta, s2, sigma2) {
    const d = 1 - (1 - eta * s2) * (1 - eta * s2);
    return d > 0 ? (eta * eta * s2 * sigma2) / d : Infinity;
  }

  /** xorshift + Box-Muller，可复现 */
  function makeRng(seed) {
    let s = (seed >>> 0) || 1, spare = null;
    function uniform() {
      s = (s + 0x6d2b79f5) >>> 0;
      let t = s;
      t = Math.imul(t ^ (t >>> 15), t | 1);
      t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    }
    function gauss() {
      if (spare !== null) { const v = spare; spare = null; return v; }
      let u = 0, v = 0;
      while (u === 0) u = uniform();
      while (v === 0) v = uniform();
      const r = Math.sqrt(-2 * Math.log(u));
      spare = r * Math.sin(2 * Math.PI * v);
      return r * Math.cos(2 * Math.PI * v);
    }
    return { uniform, gauss };
  }

  /** 返回每步的 MSD（对 paths 条独立轨道的均值），e_0 = 1 */
  function simulateLMS(etaFn, s2, sigma2, paths, steps, seed) {
    const rng = makeRng(seed === undefined ? 12345 : seed);
    const sd = Math.sqrt(sigma2);
    const e = new Float64Array(paths); e.fill(1);
    const series = new Float64Array(steps + 1);
    series[0] = 1;
    for (let k = 0; k < steps; k++) {
      const eta = etaFn(k);
      const a = 1 - eta * s2, g = eta * Math.sqrt(s2) * sd;
      let sum = 0;
      for (let i = 0; i < paths; i++) { const v = a * e[i] + g * rng.gauss(); e[i] = v; sum += v * v; }
      series[k + 1] = sum / paths;
    }
    return series;
  }

  /** 对一段序列给出均值与标准误（统计风味：给蒙特卡洛点配误差棒） */
  function windowStats(series, from, to) {
    let sum = 0, n = 0;
    for (let i = from; i < to; i++) { sum += series[i]; n++; }
    const mean = sum / n;
    let ss = 0;
    for (let i = from; i < to; i++) { const d = series[i] - mean; ss += d * d; }
    const sd = n > 1 ? Math.sqrt(ss / (n - 1)) : 0;
    return { mean, se: sd / Math.sqrt(n), n };
  }

  root.StabilityCore = {
    // 复数
    cabs, cscale, cstr,
    // 线性代数
    eig2x2, spectralRadius, isHurwitz, solveLyapunov2x2,
    // 前向 Euler
    eulerMultiplier, eulerStable, eulerRho, eulerEtaMax, bisectEtaMax, eulerEffectiveRate,
    // 模拟
    simulateDiagonal, exactFlowDiagonal, norm,
    // heavy-ball
    heavyBallMatrix, heavyBallEigs, heavyBallRho, heavyBallEtaMax, heavyBallDisc,
    continuousDamping, simulateHeavyBall,
    // 噪声
    steadyVariance, makeRng, simulateLMS, windowStats,
  };
})(typeof globalThis !== 'undefined' ? globalThis : this);
