/*!
 * test-core.js — 对 stability-core.js 的断言测试 + 数值证据表
 * 运行： node test-core.js
 * 输出的表格与 README §6、网页 §05「自检」同源同参。
 */
'use strict';
require('./stability-core.js');
const C = globalThis.StabilityCore;

let pass = 0, fail = 0;
const rows = [];

function fmt(v, d) {
  if (typeof v === 'number') {
    if (!isFinite(v)) return String(v);
    if (v !== 0 && (Math.abs(v) >= 1e4 || Math.abs(v) < 1e-3)) return v.toExponential(2);
    return v.toFixed(d === undefined ? 6 : d);
  }
  return String(v);
}
/** 断言：实测 a 与理论 b 的相对误差 < tol */
function close(name, a, b, tol, extra) {
  const err = (b === 0) ? Math.abs(a) : Math.abs(a - b) / Math.abs(b);
  const ok = err < tol;
  ok ? pass++ : fail++;
  rows.push([ok ? 'PASS' : 'FAIL', name, fmt(b), fmt(a), err.toExponential(1), extra || '']);
}
function ok(name, cond, extra) {
  cond ? pass++ : fail++;
  rows.push([cond ? 'PASS' : 'FAIL', name, '—', cond ? 'true' : 'false', '—', extra || '']);
}

/* ---------- 0. 基础：2x2 特征值 / Lyapunov ---------- */
(function () {
  const e = C.eig2x2([[0, -1], [1, 0]]);
  ok('eig2x2: 旋转矩阵 → ±i', Math.abs(e[0].im - 1) < 1e-15 && Math.abs(e[1].im + 1) < 1e-15);

  const P = C.solveLyapunov2x2([[-1, 0], [0, -4]], [[1, 0], [0, 1]]);
  ok('(5) JᵀP + PJ = −I 且 P ≻ 0', P && P.pd && P.residual < 1e-12,
    'P = diag(' + fmt(P.P[0][0], 4) + ', ' + fmt(P.P[1][1], 4) + ')');
})();

/* ---------- A. 线性增量迭代：η_max = 2/λ_max ---------- */
const LABELS = [1, 4];
const JA = LABELS.map((l) => ({ re: -l, im: 0 }));
(function () {
  const th = C.eulerEtaMax(JA).etaMax;
  const num = C.bisectEtaMax((e) => C.eulerRho(JA, e), 1e-12, 2 / 4, 60);
  close('A · η_max = 2/λ_max（解析 vs 二分）', num, th, 1e-12, 'λ = {1, 4}');
  close('A · η_max 数值', th, 0.5, 1e-15);

  const etas = [0.025, 0.075, 0.125, 0.2];
  const want = [0.975, 0.925, 0.875, 0.8];
  etas.forEach((e, i) => close('A · ρ(I−ηA), η = ' + e, C.eulerRho(JA, e), want[i], 1e-14));

  // 边界两侧：0.99η_max 收敛 / 1.01η_max 发散
  // 注意：发散时 simulateDiagonal 会在 |x| > 1e12 处提前返回，故取末元素
  const run = (eta) => { const a = C.simulateDiagonal(LABELS, eta, 2000, [1, 1]); return C.norm(a[a.length - 1]); };
  const lo = run(0.99 * th), hi = run(1.01 * th);
  ok('A · 边界两侧（0.99 收敛 / 1.01 发散）', lo < 1e-8 && hi > 1e10,
    fmt(lo, 2) + ' vs ' + fmt(hi, 2));

  // 复特征值的步长上界：2(-Reλ)/|λ|²
  const lc = { re: -1, im: 3 };
  close('A · 复 λ = −1±3i 的 η_max = 2/|λ|²', C.eulerEtaMax([lc]).etaMax, 0.2, 1e-15);
})();

/* ---------- B. Heavy-ball ---------- */
(function () {
  const lamMax = 2, beta = 0.5;
  [0, 0.2, 0.5, 0.8, 0.9].forEach((b) => {
    const th = C.heavyBallEtaMax(b, lamMax);
    const num = C.bisectEtaMax((e) => C.heavyBallRho(e, b, lamMax), 1e-12, 4, 70);
    close('B · η_max(β=' + b + ') = 2(1+β)/λ_max', num, th, 1e-12);
  });

  // 复根区 |μ| = √β，与 η 无关
  const etaC = 0.7;
  ok('B · η=' + etaC + ' 落在复根区', C.heavyBallDisc(etaC, beta, lamMax) < 0);
  close('B · 复根区 ρ = √β', C.heavyBallRho(etaC, beta, lamMax), Math.sqrt(beta), 1e-14);

  // 实测几何增长率 vs 谱半径
  const traj = C.simulateHeavyBall(etaC, beta, lamMax, 1800, 1);
  const n0 = Math.abs(traj[600]), n1 = Math.abs(traj[1800]);
  const meas = Math.exp((Math.log(n1) - Math.log(n0)) / 1200);
  close('B · 实测增长率 vs ρ（几何平均）', meas, C.heavyBallRho(etaC, beta, lamMax), 1e-3);

  // GAP：连续稳定 + 离散爆炸
  const etaGap = 1.5 * C.heavyBallEtaMax(beta, lamMax);
  const a = C.continuousDamping(etaGap, beta);
  const t2 = C.simulateHeavyBall(etaGap, beta, lamMax, 200, 1);
  const blow = t2.findIndex((v) => Math.abs(v) > 1e12);
  ok('B · GAP：ODE 阻尼 a > 0 而离散爆炸', a > 0 && blow > 0 && blow < 60,
    'a=' + fmt(a, 4) + '，第 ' + blow + ' 步 |x|>1e12，ρ=' + fmt(C.heavyBallRho(etaGap, beta, lamMax), 4));
})();

/* ---------- C. 噪声地板 ---------- */
(function () {
  const eta = 0.05, s2 = 1, sigma2 = 1;
  const th = C.steadyVariance(eta, s2, sigma2);
  const s = C.simulateLMS(() => eta, s2, sigma2, 4000, 1500, 20240607);
  const w = C.windowStats(s, 1000, 1500);
  close('C · 稳态方差 V = ησ²/(2−ηs²) vs MC', w.mean, th, 0.02,
    '±1.96SE = ' + fmt(1.96 * w.se, 5) + '（时间序列自相关，SE 偏小）');

  // 小 η 渐近：V ≈ ησ²/2
  const e2 = 0.002;
  close('C · 小 η 渐近 V/(ησ²/2) → 1', C.steadyVariance(e2, s2, sigma2) / (e2 * sigma2 / 2), 1, 2e-3);

  // 蒙特卡洛对数斜率 ≈ 1（MSD ∝ η）
  const pts = [0.01, 0.02, 0.04, 0.08].map((e, i) => {
    const ser = C.simulateLMS(() => e, s2, sigma2, 2000, 900, 500 + i);
    return [Math.log(e), Math.log(C.windowStats(ser, 600, 900).mean)];
  });
  const slope = (pts[3][1] - pts[0][1]) / (pts[3][0] - pts[0][0]);
  close('C · log-log 斜率 d log MSD / d log η', slope, 1, 0.03);

  // 衰减步长 vs 常数步长：η_k = c/(k+10) ⟹ V ≈ c/(2k)，按 1/k 下降
  const decay = C.simulateLMS((k) => 2 / (k + 10), s2, sigma2, 2000, 3000, 909);
  const floor = C.steadyVariance(0.2, s2, sigma2);
  ok('C · 衰减步长远低于同期的常步长噪声地板', decay[3000] < 0.01 * floor,
    'MSD(3000) = ' + fmt(decay[3000], 2) + '，地板 = ' + fmt(floor, 4));
  close('C · 衰减步长 MSD·k → c/2（1/k 律）', decay[3000] * 3000, 1, 0.35);

  // 前向 Euler 全局误差 O(η)：err/η 收敛到一个常数
  const euler = (eta2, T) => { let x = 1; const n = Math.round(T / eta2); for (let i = 0; i < n; i++) x *= (1 - eta2); return Math.abs(x - Math.exp(-T)); };
  const ratios = [0.05, 0.01, 0.002].map((e) => euler(e, 1) / e);
  close('C · 前向 Euler 全局误差 = O(η)', ratios[2], ratios[1], 0.01,
    'err/η = ' + ratios.map((v) => fmt(v, 5)).join(' / '));
})();

/* ---------- 输出 ---------- */
const W = [6, 46, 18, 18, 10];
const line = (a) => a[0].padEnd(W[0]) + a[1].padEnd(W[1]) + a[2].padStart(W[2]) + a[3].padStart(W[3]) + a[4].padStart(W[4]) + (a[5] ? '   ' + a[5] : '');
console.log('stability-core.js 断言与数值证据');
console.log('='.repeat(118));
console.log(line(['状态', '检验项', '理论值', '实测值', '相对误差', '备注']));
console.log('-'.repeat(118));
rows.forEach((r) => console.log(line(r)));
console.log('='.repeat(118));
console.log(pass + ' passed, ' + fail + ' failed, ' + rows.length + ' checks');
process.exit(fail ? 1 : 0);
