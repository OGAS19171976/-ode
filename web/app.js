/*!
 * app.js — 交互式教学页的渲染与交互层（依赖 stability-core.js）
 * 四个演示：
 *   1. 时间重标定：离散迭代点是否落在 ODE 流上
 *   2. 稳定域：左半平面 vs 以 -1/eta 为圆心的圆盘（拖拽 lambda 点）
 *   3. Hopf 分岔：Heavy-ball 特征值穿越单位圆
 *   4. 噪声地板：稳态 MSD ∝ eta
 */
'use strict';
(function () {
  const C = globalThis.StabilityCore;

  const T = {
    bg: '#0e1523', panel: '#0b111d', grid: '#1a2334', axis: '#38455f',
    text: '#c9d4e8', dim: '#7b8aa5', accent: '#4da3ff', good: '#3ddc97',
    bad: '#ff6b6b', warn: '#ffb454', curve: '#8b7bff', soft: '#5a6b88',
  };

  /* ============================ 通用绘图 ============================ */
  function niceTicks(min, max, n) {
    if (!(max > min)) return [min];
    const span = max - min, raw = span / Math.max(1, n);
    const mag = Math.pow(10, Math.floor(Math.log10(raw))), norm = raw / mag;
    const step = mag * (norm < 1.5 ? 1 : norm < 3 ? 2 : norm < 7 ? 5 : 10);
    const out = [];
    for (let t = Math.ceil(min / step) * step; t <= max + span * 1e-9; t += step) out.push(Math.abs(t) < 1e-12 ? 0 : t);
    return out;
  }
  function fmtTick(v) {
    const a = Math.abs(v);
    if (a === 0) return '0';
    if (a >= 1e4 || a < 1e-3) return v.toExponential(0);
    return String(Math.round(v * 1e6) / 1e6);
  }
  const pow10fmt = (v) => '1e' + Math.round(v);

  class Plot {
    constructor(canvas) {
      this.cv = canvas; this.ctx = canvas.getContext('2d');
      this.pad = { l: 54, r: 16, t: 16, b: 40 };
      this.world = { xmin: 0, xmax: 1, ymin: 0, ymax: 1 };
      this.W = 1; this.H = 1;
    }
    get box() {
      return {
        x: this.pad.l, y: this.pad.t,
        w: Math.max(1, this.W - this.pad.l - this.pad.r),
        h: Math.max(1, this.H - this.pad.t - this.pad.b),
      };
    }
    resize() {
      const dpr = Math.min(window.devicePixelRatio || 1, 2);
      const w = Math.max(1, this.cv.clientWidth), h = Math.max(1, this.cv.clientHeight);
      if (this.cv.width !== Math.round(w * dpr) || this.cv.height !== Math.round(h * dpr)) {
        this.cv.width = Math.round(w * dpr); this.cv.height = Math.round(h * dpr);
      }
      this.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      this.W = w; this.H = h;
      return this;
    }
    setWorld(xmin, xmax, ymin, ymax) { this.world = { xmin, xmax, ymin, ymax }; return this; }
    /** y 范围固定，x 范围由画布宽高比推出 → 保证圆是圆 */
    setEqualAspect(ymin, ymax, anchorX, anchorFrac) {
      const b = this.box, xr = (ymax - ymin) * (b.w / b.h);
      const xmax = anchorX + anchorFrac * xr;
      this.world = { xmin: xmax - xr, xmax, ymin, ymax };
      return this;
    }
    X(x) { const b = this.box, w = this.world; return b.x + (x - w.xmin) / (w.xmax - w.xmin) * b.w; }
    Y(y) { const b = this.box, w = this.world; return b.y + b.h - (y - w.ymin) / (w.ymax - w.ymin) * b.h; }
    invX(px) { const b = this.box, w = this.world; return w.xmin + (px - b.x) / b.w * (w.xmax - w.xmin); }
    invY(py) { const b = this.box, w = this.world; return w.ymin + (b.y + b.h - py) / b.h * (w.ymax - w.ymin); }
    scaleX() { return this.box.w / (this.world.xmax - this.world.xmin); }

    drawFrame(o) {
      o = o || {};
      const c = this.ctx, b = this.box, w = this.world;
      const xt = niceTicks(w.xmin, w.xmax, o.nx || 5);
      const yt = niceTicks(w.ymin, w.ymax, o.ny || 5);
      const xf = o.xfmt || fmtTick, yf = o.yfmt || fmtTick;
      c.save();
      c.fillStyle = T.bg; c.fillRect(0, 0, this.W, this.H);
      c.fillStyle = T.panel; c.fillRect(b.x, b.y, b.w, b.h);
      c.strokeStyle = T.grid; c.lineWidth = 1;
      for (const t of xt) { const px = Math.round(this.X(t)) + .5; c.beginPath(); c.moveTo(px, b.y); c.lineTo(px, b.y + b.h); c.stroke(); }
      for (const t of yt) { const py = Math.round(this.Y(t)) + .5; c.beginPath(); c.moveTo(b.x, py); c.lineTo(b.x + b.w, py); c.stroke(); }
      if (o.zeroAxes !== false) {
        c.strokeStyle = T.axis; c.lineWidth = 1.2;
        if (w.ymin < 0 && w.ymax > 0) { const py = Math.round(this.Y(0)) + .5; c.beginPath(); c.moveTo(b.x, py); c.lineTo(b.x + b.w, py); c.stroke(); }
        if (w.xmin < 0 && w.xmax > 0) { const px = Math.round(this.X(0)) + .5; c.beginPath(); c.moveTo(px, b.y); c.lineTo(px, b.y + b.h); c.stroke(); }
      }
      c.font = '11px ui-monospace, Menlo, Consolas, monospace';
      c.fillStyle = T.dim; c.textAlign = 'center'; c.textBaseline = 'top';
      for (const t of xt) c.fillText(xf(t), this.X(t), b.y + b.h + 6);
      c.textAlign = 'right'; c.textBaseline = 'middle';
      for (const t of yt) c.fillText(yf(t), b.x - 8, this.Y(t));
      c.fillStyle = T.text; c.font = '12px system-ui, sans-serif';
      if (o.xlabel) { c.textAlign = 'center'; c.textBaseline = 'bottom'; c.fillText(o.xlabel, b.x + b.w / 2, this.H - 2); }
      if (o.ylabel) {
        c.save(); c.translate(13, b.y + b.h / 2); c.rotate(-Math.PI / 2);
        c.textAlign = 'center'; c.textBaseline = 'top'; c.fillText(o.ylabel, 0, 0); c.restore();
      }
      c.restore();
      return this;
    }
    clip() { const c = this.ctx, b = this.box; c.save(); c.beginPath(); c.rect(b.x, b.y, b.w, b.h); c.clip(); }
    unclip() { this.ctx.restore(); }

    polyline(pts, o) {
      o = o || {};
      const c = this.ctx;
      c.save(); c.strokeStyle = o.color || T.accent; c.lineWidth = o.width || 2;
      if (o.dash) c.setLineDash(o.dash);
      c.beginPath();
      let started = false;
      for (const p of pts) {
        if (!isFinite(p[0]) || !isFinite(p[1])) { started = false; continue; }
        const x = this.X(p[0]), y = this.Y(p[1]);
        if (!started) { c.moveTo(x, y); started = true; } else c.lineTo(x, y);
      }
      c.stroke(); c.restore();
      return this;
    }
    points(pts, o) {
      o = o || {};
      const c = this.ctx;
      c.save(); c.fillStyle = o.color || T.good;
      for (const p of pts) {
        if (!isFinite(p[0]) || !isFinite(p[1])) continue;
        const x = this.X(p[0]), y = this.Y(p[1]);
        if (x < this.box.x - 40 || x > this.box.x + this.box.w + 40 || y < this.box.y - 40 || y > this.box.y + this.box.h + 40) continue;
        c.beginPath(); c.arc(x, y, o.r || 2.6, 0, 6.2832); c.fill();
      }
      c.restore();
      return this;
    }
    dot(x, y, o) {
      o = o || {};
      const c = this.ctx, px = this.X(x), py = this.Y(y);
      c.save();
      c.beginPath(); c.arc(px, py, o.r || 4, 0, 6.2832);
      c.fillStyle = o.color || '#fff'; c.fill();
      if (o.ring) { c.strokeStyle = o.ring; c.lineWidth = 2; c.stroke(); }
      c.restore();
      return this;
    }
    hline(y, o) { o = o || {}; const b = this.box; return this.polyline([[this.world.xmin, y], [this.world.xmax, y]], o); }
    vline(x, o) { o = o || {}; return this.polyline([[x, this.world.ymin], [x, this.world.ymax]], o); }
    circle(cx, cy, r, o) {
      o = o || {};
      const c = this.ctx, s = this.scaleX(), b = this.box;
      c.save();
      c.beginPath(); c.arc(this.X(cx), this.Y(cy), r * s, 0, 6.2832);
      if (o.fill) { c.fillStyle = o.fill; c.fill(); }
      if (o.stroke) { c.strokeStyle = o.stroke; c.lineWidth = o.width || 1.6; if (o.dash) c.setLineDash(o.dash); c.stroke(); }
      c.restore();
      // 标注
      if (o.label) {
        const px = this.X(cx), py = this.Y(cy);
        c.save(); c.fillStyle = T.dim; c.font = '11px system-ui, sans-serif';
        c.textAlign = 'center'; c.textBaseline = 'top';
        c.fillText(o.label, Math.min(Math.max(px, b.x + 40), b.x + b.w - 40), Math.min(py + r * s + 4, b.y + b.h - 14));
        c.restore();
      }
      return this;
    }
    text(x, y, s, o) {
      o = o || {};
      const c = this.ctx;
      c.save(); c.fillStyle = o.color || T.text;
      c.font = o.font || '12px system-ui, sans-serif';
      c.textAlign = o.align || 'left'; c.textBaseline = o.baseline || 'alphabetic';
      c.fillText(s, this.X(x), this.Y(y));
      c.restore();
      return this;
    }
    label(px, py, s, o) {
      o = o || {};
      const b = this.box, c = this.ctx;
      c.save(); c.fillStyle = o.color || T.text;
      c.font = o.font || '12px system-ui, sans-serif';
      c.textAlign = o.align || 'left'; c.textBaseline = o.baseline || 'alphabetic';
      c.fillText(s, Math.min(Math.max(px, b.x + 2), b.x + b.w - 2), Math.min(Math.max(py, b.y + 12), b.y + b.h - 2));
      c.restore();
      return this;
    }
    errorBar(x, y, half, o) {
      o = o || {};
      const c = this.ctx, px = this.X(x);
      c.save(); c.strokeStyle = o.color || T.good; c.lineWidth = 1.2;
      c.beginPath();
      c.moveTo(px, this.Y(y - half)); c.lineTo(px, this.Y(y + half));
      c.moveTo(px - 4, this.Y(y - half)); c.lineTo(px + 4, this.Y(y - half));
      c.moveTo(px - 4, this.Y(y + half)); c.lineTo(px + 4, this.Y(y + half));
      c.stroke(); c.restore();
      return this;
    }
    legend(items, o) {
      o = o || {};
      const c = this.ctx, b = this.box;
      const x = o.left ? b.x + 12 : b.x + b.w - 12;
      const align = o.left ? 'left' : 'right';
      c.save(); c.font = '11.5px system-ui, sans-serif'; c.textAlign = align; c.textBaseline = 'middle';
      items.forEach((it, i) => {
        const y = b.y + 14 + i * 17;
        const tx = o.left ? x + 18 : x - 18;
        c.fillStyle = it.color;
        if (it.dot) { c.beginPath(); c.arc(o.left ? x + 5 : x - 5, y, 3.4, 0, 6.2832); c.fill(); }
        else { c.fillRect(o.left ? x : x - 12, y - 1.6, 12, 3.2); }
        c.fillStyle = T.dim; c.fillText(it.label, tx, y);
      });
      c.restore();
      return this;
    }
    clearCanvas() { this.ctx.save(); this.ctx.setTransform(1, 0, 0, 1, 0, 0); this.ctx.clearRect(0, 0, this.cv.width, this.cv.height); this.ctx.restore(); return this; }
  }

  const $ = (s) => document.querySelector(s);
  const $$ = (s) => Array.prototype.slice.call(document.querySelectorAll(s));
  const fx = (v, d) => Number(v).toFixed(d === undefined ? 3 : d);
  const plots = {};

  /* ============================ 演示 1 ============================ */
  const D1 = { eta: 0.15, lams: [1, 4], rescale: true, K: 40 };

  function renderD1() {
    const eta = D1.eta, lams = D1.lams, K = D1.K;
    const J = lams.map((l) => ({ re: -l, im: 0 }));
    const info = C.eulerEtaMax(J);
    const etaMax = info.etaMax;
    const rho = Math.max.apply(null, lams.map((l) => Math.abs(1 - eta * l)));
    const stable = rho < 1;
    const xk = C.simulateDiagonal(lams, eta, K, [1, 1]);
    const tend = K * eta;
    const diverge = Math.abs(xk[xk.length - 1][1]) > 1e6;

    const a = plots.p1a.resize();
    a.setWorld(-1.05, 1.15, -2.2, 1.3);
    a.drawFrame({ xlabel: 'x₁', ylabel: 'x₂', nx: 4, ny: 5 });
    a.clip();
    const cont = [];
    for (let i = 0; i <= 320; i++) cont.push(C.exactFlowDiagonal(lams, tend * i / 320, [1, 1]));
    a.polyline(cont, { color: T.accent, width: 2.6 });
    a.points(xk.map((p) => [p[0], p[1]]), { color: stable ? T.good : T.bad, r: 2.7 });
    a.dot(1, 1, { color: '#ffffff', r: 3.6 });
    a.dot(0, 0, { color: T.warn, r: 3.2 });
    a.unclip();
    a.legend([
      { label: 'ODE 精确流 x(t)', color: T.accent },
      { label: '离散迭代 x_k（同一张图上）', color: stable ? T.good : T.bad, dot: true },
    ], { left: true });

    const b = plots.p1b.resize();
    b.setWorld(0, K, -20, 3.2);
    b.drawFrame({ xlabel: '迭代步 k', ylabel: 'log₁₀‖x_k − x*‖', yfmt: pow10fmt, nx: 5, ny: 5 });
    b.clip();
    const contPts = [];
    for (let k = 0; k <= K; k++) {
      const t = D1.rescale ? k * eta : k;
      const n = C.norm(C.exactFlowDiagonal(lams, t, [1, 1]));
      contPts.push([k, n > 0 ? Math.log10(n) : -20]);
    }
    b.polyline(contPts, { color: T.accent, width: 2.4, dash: D1.rescale ? null : [6, 4] });
    const dp = [];
    xk.forEach((x, k) => { const n = C.norm(x); if (k <= K) dp.push([k, n > 0 ? Math.max(Math.log10(n), -20) : -20]); });
    b.points(dp, { color: stable ? T.good : T.bad, r: 2.4 });
    b.unclip();
    b.legend([
      { label: D1.rescale ? 'ODE 流 @ t = kη（重标定）' : 'ODE 流 @ t = k（未重标定）', color: T.accent },
      { label: '离散迭代', color: stable ? T.good : T.bad, dot: true },
    ], { left: true });

    const rateC = lams[0];
    const rateD = C.eulerEffectiveRate({ re: -lams[0], im: 0 }, eta);
    $('#d1-info').innerHTML = [
      ['η', fx(eta, 3) + '   (η_max = ' + fx(etaMax, 4) + ')'],
      ['谱半径 ρ(I − ηA) = maxᵢ|1 − ηλᵢ|', fx(rho, 6)],
      ['状态', stable ? '<span class="tag ok">离散收敛</span>' : '<span class="tag bad">离散发散</span>'],
      ['最慢模态（λ_min = 1）连续速率', fx(rateC, 6)],
      ['同一模态的离散有效速率 −ln ρ/η', fx(rateD, 6)],
      ['离散相对延迟 ≈ ηλ/2', (100 * (rateC - rateD) / rateC).toFixed(2) + '%'],
    ].map((r) => '<tr><td>' + r[0] + '</td><td>' + r[1] + '</td></tr>').join('');

    $('#d1-tip').innerHTML = D1.rescale
      ? (diverge || !stable
        ? 'η 已越过 η_max = 2/λ_max：连续时间不可能不稳定，但 <b>离散格式的绝对稳定域是有限的</b>，迭代开始振荡放大。'
        : '重标定后 t = kη，离散点整条贴在 ODE 流上；η 越小越贴 —— 这正是"离散算法 = ODE 的一个数值积分格式"。')
      : '关掉重标定：把第 k 步画在 t = k。两条曲线立刻相差一个 <b>η</b> 的时间尺度 —— 轨迹并没有错，只是时间轴没换算。这是最容易误判"ODE 与仿真不一致"的原因。';
  }

  /* ============================ 演示 2 ============================ */
  const D2 = {
    eta: 0.30,
    pts: [
      { re: -0.5, im: 0 }, { re: -2, im: 0 }, { re: -4, im: 0 },
      { re: -1, im: 3 }, { re: -0.3, im: 3.6 }, { re: 1.5, im: 0 },
    ],
    drag: -1,
  };
  const D2_PRESETS = {
    real: [{ re: 1, im: 0 }, { re: 2, im: 0 }, { re: 4, im: 0 }, { re: 0.5, im: 0 }],
    complex: [{ re: -1, im: 3 }, { re: -0.3, im: 3.6 }, { re: -2, im: 1.2 }, { re: -0.2, im: 1 }],
    ill: [{ re: -0.05, im: 0 }, { re: -4, im: 0 }, { re: -0.2, im: 2.5 }, { re: -3, im: 1 }],
  };

  function d2Verdict(l, eta) {
    const m = Math.hypot(1 + eta * l.re, eta * l.im);
    if (m < 1) return { cls: 'ok', text: '离散稳定' };
    if (l.re < 0) return { cls: 'warn', text: '连续稳 / 离散不稳（缝隙）' };
    return { cls: 'bad', text: '连续与离散都不稳' };
  }

  function renderD2() {
    const eta = D2.eta, p = plots.p2.resize();
    p.setEqualAspect(-5.6, 5.6, 0, 0.20);
    p.drawFrame({ xlabel: 'Re λ', ylabel: 'Im λ', nx: 8, ny: 6 });
    p.clip();
    // 连续时间稳定区：左半平面
    p.ctx.save();
    p.ctx.fillStyle = 'rgba(61,220,151,.055)';
    p.ctx.fillRect(p.X(p.world.xmin), p.Y(0), p.X(0) - p.X(p.world.xmin), p.Y(p.world.ymin) - p.Y(0));
    p.ctx.restore();
    // 前向 Euler 的绝对稳定域：|1+eta*lam|<1 ⇔ |lam + 1/eta| < 1/eta
    p.circle(-1 / eta, 0, 1 / eta, { fill: 'rgba(77,163,255,.10)', stroke: T.accent, width: 1.6, label: '|1+ηλ| < 1' });
    p.text(p.world.xmin + 0.4, 4.6, '左半平面 = 连续时间稳定', { color: 'rgba(61,220,151,.85)' });
    // λ 点
    D2.pts.forEach((l, i) => {
      const v = d2Verdict(l, eta);
      const col = v.cls === 'ok' ? T.good : v.cls === 'warn' ? T.warn : T.bad;
      const r = 4.4 * (D2.drag === i ? 1.35 : 1);
      p.dot(l.re, l.im, { color: col, r });
      p.label(p.X(l.re) + 9, p.Y(l.im) - 8, 'λ' + (i + 1), { color: T.dim, font: '11px ui-monospace, monospace' });
    });
    p.unclip();

    $('#d2-table').innerHTML = '<table class="readout">' + D2.pts.map((l, i) => {
      const m = Math.hypot(1 + eta * l.re, eta * l.im);
      const v = d2Verdict(l, eta);
      return '<tr><td>λ' + (i + 1) + ' = ' + fx(l.re, 2) + (l.im >= 0 ? ' + ' : ' − ') + fx(Math.abs(l.im), 2) + 'i</td><td>' +
        fx(m, 4) + ' <span class="tag ' + v.cls + '">' + v.text + '</span></td></tr>';
    }).join('') + '</table>';

    const gap = D2.pts.filter((l) => l.re < 0 && Math.hypot(1 + eta * l.re, eta * l.im) >= 1);
    $('#d2-tip').innerHTML = '稳定判据不是一个半平面，而是一个<b>圆盘</b>：<code>|1+ηλ| &lt; 1 ⟺ |λ + 1/η| &lt; 1/η</code>。' +
      '拖动任意 λ 点试试。当前 η = ' + fx(eta, 3) + '，圆盘圆心 ' + fx(-1 / eta, 2) + '、半径 ' + fx(1 / eta, 2) +
      (gap.length
        ? '。<b>' + gap.length + ' 个点落在缝隙里</b>：它们的 ODE 是稳定的，离散迭代却是发散的 —— 调小学习率会放大圆盘把它们救回来，但特征值符号（结构性不稳定）救不回来。'
        : '。当前所有 λ 都在圆盘内：连续稳定 ⟹ 离散稳定。把 η 调大，圆盘收缩，缝隙就出现了。');
  }

  /* ============================ 演示 3 ============================ */
  const D3 = { eta: 0.6, beta: 0.5, lam: 2 };

  function renderD3() {
    const eta = D3.eta, beta = D3.beta, lam = D3.lam;
    const eigs = C.heavyBallEigs(eta, beta, lam);
    const rho = C.spectralRadius(eigs);
    const disc = C.heavyBallDisc(eta, beta, lam);
    const isComplex = disc < 0;
    const etaMax = C.heavyBallEtaMax(beta, lam);
    const a = C.continuousDamping(eta, beta);
    const stable = rho < 1;

    const p = plots.p3a.resize();
    p.setEqualAspect(-2.7, 2.7, 0, 0.22);
    p.drawFrame({ xlabel: 'Re μ', ylabel: 'Im μ', nx: 6, ny: 5 });
    p.clip();
    p.circle(0, 0, 1, { stroke: T.soft, width: 1.3, dash: [5, 4], label: '单位圆 |μ| = 1' });
    // η 扫描的根轨迹
    const locus = [];
    for (let e = 0.02; e <= 3.2; e += 0.012) locus.push(C.heavyBallEigs(e, beta, lam));
    locus.forEach((pair) => { p.polyline(pair.map((z) => [z.re, z.im]), { color: 'rgba(139,123,255,.35)', width: 1.2 }); });
    eigs.forEach((z) => {
      p.dot(z.re, z.im, { color: stable ? T.good : T.bad, r: 4.6, ring: 'rgba(255,255,255,.25)' });
    });
    // 由根到单位圆的连线，直观读出模长
    eigs.forEach((z) => p.polyline([[0, 0], [z.re, z.im]], { color: 'rgba(255,255,255,.16)', width: 1 }));
    p.unclip();

    const q = plots.p3b.resize();
    q.setWorld(0, 60, -18, 11);
    q.drawFrame({ xlabel: '迭代步 k', ylabel: 'log₁₀|x_k|', yfmt: pow10fmt, nx: 6, ny: 6 });
    q.clip();
    const traj = C.simulateHeavyBall(eta, beta, lam, 60, 1);
    const tp = [];
    traj.forEach((v, k) => { if (k <= 60) tp.push([k, Math.abs(v) > 0 ? Math.max(Math.log10(Math.abs(v)), -18) : -18]); });
    q.polyline(tp, { color: stable ? T.good : T.bad, width: 2.2 });
    if (stable) q.hline(Math.log10(1e-16), { color: T.grid, width: 1 });
    q.unclip();
    q.legend([{ label: stable ? (isComplex ? '振荡收敛（复特征值）' : '单调收敛') : '发散', color: stable ? T.good : T.bad }], { left: true });

    const status = !stable
      ? '<span class="tag bad">' + (isComplex ? '振荡发散' : '单调发散') + '</span>'
      : (isComplex ? '<span class="tag ok">振荡收敛</span>' : '<span class="tag ok">单调收敛</span>');
    $('#d3-info').innerHTML = [
      ['β', fx(beta, 2)],
      ['η', fx(eta, 3) + '   (η_max = 2(1+β)/λ_max = ' + fx(etaMax, 4) + ')'],
      ['特征值 μ₁,₂', C.cstr(eigs[0], 3) + ' , ' + C.cstr(eigs[1], 3)],
      ['判别式 tr²−4det', fx(disc, 4) + (isComplex ? ' < 0（复根，|μ| = √β）' : ' ≥ 0（实根）')],
      ['谱半径 ρ = max|μ|', fx(rho, 6)],
      ['状态', status],
      ['ODE 阻尼 a = (1−β)/√η', fx(a, 4) + (a > 0 ? ' > 0 ⟹ ODE 无条件稳定' : ' ≤ 0 ⟹ ODE 本身不稳')],
    ].map((r) => '<tr><td>' + r[0] + '</td><td>' + r[1] + '</td></tr>').join('');

    $('#d3-tip').innerHTML = (a > 0 && !stable)
      ? '<b>缝隙现场：</b>连续时间阻尼 a = ' + fx(a, 3) + ' &gt; 0，ODE 是稳定的；但 η = ' + fx(eta, 2) + ' &gt; η_max = ' + fx(etaMax, 3) +
        '，离散迭代已经越过单位圆。这就是"连续稳定 ⟹ 离散稳定"不成立的构造性反例。'
      : (isComplex
        ? '复特征值区：两根共轭，模长恒为 <b>√β</b>（与 η 无关），迭代以固定几何速率振荡收敛。动量越大，β 越接近 1，振荡越慢衰减。'
        : '实特征值区：按 <code>0 &lt; ηλ &lt; 2(1+β)</code> 判据即可，上界随 β 线性放大 —— 动量用一个"更慢的振荡"换来了更大的可用步长。');
  }

  /* ============================ 演示 4 ============================ */
  const D4 = { eta: 0.15, s2: 1, sigma2: 1, paths: 800, cache: {} };

  function d4Sweep() {
    if (D4.cache.sweep) return D4.cache.sweep;
    const out = [];
    for (let i = 0; i <= 9; i++) {
      const eta = Math.pow(10, -2.3 + i * (-0.35 + 2.3) / 9);
      const s = C.simulateLMS(() => eta, D4.s2, D4.sigma2, 600, 400, 1000 + i);
      const st = C.windowStats(s, 250, 400);
      out.push({ eta, mean: st.mean, se: st.se, theory: C.steadyVariance(eta, D4.s2, D4.sigma2) });
    }
    D4.cache.sweep = out;
    return out;
  }
  function d4Series(eta) {
    const key = eta.toFixed(4);
    if (D4.cache.series && D4.cache.series.key === key) return D4.cache.series;
    const constSeries = C.simulateLMS(() => eta, D4.s2, D4.sigma2, D4.paths, 1200, 777);
    const decaySeries = C.simulateLMS((k) => 2 / (k + 10), D4.s2, D4.sigma2, D4.paths, 1200, 778);
    D4.cache.series = { key, constSeries, decaySeries };
    return D4.cache.series;
  }

  function renderD4() {
    const eta = D4.eta, s2 = D4.s2, sigma2 = D4.sigma2;
    const sweep = d4Sweep();

    const p = plots.p4a.resize();
    p.setWorld(-2.4, -0.35, -5, 0.2);
    p.drawFrame({ xlabel: 'log₁₀ η', ylabel: 'log₁₀ 稳态 MSD', xfmt: pow10fmt, yfmt: pow10fmt, nx: 5, ny: 5 });
    p.clip();
    const line = [];
    for (let i = 0; i <= 60; i++) {
      const e = Math.pow(10, -2.4 + i * (2.05) / 60);
      const v = C.steadyVariance(e, s2, sigma2);
      line.push([Math.log10(e), Math.log10(v)]);
    }
    p.polyline(line, { color: T.accent, width: 2.2 });
    sweep.forEach((d) => {
      const x = Math.log10(d.eta), y = Math.log10(d.mean);
      const lo = Math.max(d.mean - 1.96 * d.se, 1e-12), hi = d.mean + 1.96 * d.se;
      p.errorBar(x, y, (Math.log10(hi) - Math.log10(lo)) / 2, { color: T.good });
      p.dot(x, y, { color: T.good, r: 3.2 });
    });
    // 当前 eta 的位置
    const vNow = C.steadyVariance(eta, s2, sigma2);
    p.dot(Math.log10(eta), Math.log10(vNow), { color: T.warn, r: 4.5 });
    p.unclip();
    p.legend([
      { label: '理论 ησ²/(2−ηs²)', color: T.accent },
      { label: '蒙特卡洛（±1.96 SE）', color: T.good, dot: true },
    ], { left: true });

    const ser = d4Series(eta);
    const q = plots.p4b.resize();
    q.setWorld(0, 1200, -7, 0.4);
    q.drawFrame({ xlabel: '迭代步 k', ylabel: 'log₁₀ MSD', yfmt: pow10fmt, nx: 6, ny: 5 });
    q.clip();
    const mk = (s) => { const out = []; for (let k = 0; k <= 1200; k += 2) out.push([k, Math.log10(Math.max(s[k], 1e-12))]); return out; };
    q.polyline(mk(ser.constSeries), { color: T.warn, width: 2 });
    q.polyline(mk(ser.decaySeries), { color: T.good, width: 2 });
    q.hline(Math.log10(vNow), { color: 'rgba(255,180,84,.5)', width: 1.2, dash: [6, 4] });
    q.unclip();
    q.legend([
      { label: '常数步长 η = ' + fx(eta, 3), color: T.warn },
      { label: '衰减步长 η_k = 2/(k+10)', color: T.good },
    ], { left: true });

    const st = C.windowStats(ser.constSeries, 900, 1200);
    $('#d4-info').innerHTML = [
      ['常步长理论稳态方差', fx(vNow, 5)],
      ['常步长蒙特卡洛（末 300 步）', fx(st.mean, 5) + ' ± ' + fx(1.96 * st.se, 5)],
      ['相对误差', (100 * Math.abs(st.mean - vNow) / vNow).toFixed(2) + '%'],
      ['小 η 渐近 ησ²/2', fx(eta * sigma2 / 2, 5)],
      ['衰减步长末值 MSD', ser.decaySeries[1200].toExponential(2)],
      ['Ση = ∞, Ση² < ∞', '几乎必然收敛到点（无噪声地板）'],
    ].map((r) => '<tr><td>' + r[0] + '</td><td>' + r[1] + '</td></tr>').join('');

    $('#d4-tip').innerHTML = '噪声<b>不移动</b>零点（h(x*) = 0 与噪声无关），它只决定误差停在哪里：常数步长的稳态 MSD 正比于 η（双对数图上斜率 1），于是误差永久停在 <code>O(√η)</code> 的噪声地板上；换成衰减步长就从"停住"变成"收敛到点"。';
  }

  /* ============================ 页内自检 ============================ */
  function renderSelfTest() {
    const rows = [];
    const lams = [1, 4];
    const J = lams.map((l) => ({ re: -l, im: 0 }));
    const th = C.eulerEtaMax(J).etaMax;
    const num = C.bisectEtaMax((e) => C.eulerRho(J, e), 1e-12, 2 / Math.max.apply(null, lams), 60);
    rows.push(['A · η_max = 2/λ_max', fx(th, 8), fx(num, 8), (Math.abs(num - th) / th).toExponential(1)]);

    const rhoVals = [0.025, 0.075, 0.125, 0.2];
    const thRho = rhoVals.map((e) => C.eulerRho(J, e));
    const wantRho = [0.975, 0.925, 0.875, 0.8];
    rows.push(['A · ρ(I−ηA) 四个 η', wantRho.join(' / '), thRho.map((v) => fx(v, 8)).join(' / '),
      Math.max.apply(null, thRho.map((v, i) => Math.abs(v - wantRho[i]))).toExponential(1)]);

    const b = 0.5, lamB = 2;
    const thB = C.heavyBallEtaMax(b, lamB);
    const numB = C.bisectEtaMax((e) => C.heavyBallRho(e, b, lamB), 1e-12, 4, 70);
    rows.push(['B · η_max(β=0.5), λ_max=2', fx(thB, 8), fx(numB, 8), (Math.abs(numB - thB) / thB).toExponential(1)]);

    const eC = 1.0;
    const rhoC = C.heavyBallRho(eC, b, lamB);
    rows.push(['B · 复根区 ρ = √β', fx(Math.sqrt(b), 8), fx(rhoC, 8), Math.abs(rhoC - Math.sqrt(b)).toExponential(1)]);

    const etaGap = 1.5 * thB, aGap = C.continuousDamping(etaGap, b);
    const traj = C.simulateHeavyBall(etaGap, b, lamB, 200, 1);
    const blow = traj.findIndex((v) => Math.abs(v) > 1e12);
    rows.push(['B · GAP: ODE 阻尼 a > 0', fx(aGap, 4) + ' > 0', blow > 0 ? '第 ' + blow + ' 步 |x| > 1e12' : '未爆', '—']);

    const P = C.solveLyapunov2x2([[-1, 0], [0, -4]], [[1, 0], [0, 1]]);
    rows.push(['(5) Lyapunov JᵀP+PJ = −I', 'P ≻ 0', P ? (P.pd ? 'P 正定' : 'P 非正定') + '，残差 ' + P.residual.toExponential(1) : '—', '—']);

    const etaN = 0.05, vTh = C.steadyVariance(etaN, 1, 1);
    const s = C.simulateLMS(() => etaN, 1, 1, 2000, 800, 20240607);
    const w = C.windowStats(s, 500, 800);
    rows.push(['C · 稳态方差 理论 vs MC', fx(vTh, 5), fx(w.mean, 5), (Math.abs(w.mean - vTh) / vTh).toExponential(1)]);

    const euler = (eta, T0) => { let x = 1; const n = Math.round(T0 / eta); for (let i = 0; i < n; i++) x *= (1 - eta); return Math.abs(x - Math.exp(-T0)); };
    const ratios = [0.05, 0.01, 0.002].map((e) => euler(e, 1) / e);
    rows.push(['C · 前向 Euler 全局误差 O(η)', 'err/η → const', ratios.map((v) => fx(v, 5)).join(' / '),
      (Math.abs(ratios[0] - ratios[2]) / ratios[2]).toExponential(1)]);

    $('#selftest').innerHTML = '<table class="readout" style="font-size:13px">' +
      '<tr><td style="color:#7b8aa5">检验项</td><td style="text-align:left;color:#7b8aa5">理论值</td>' +
      '<td style="text-align:left;color:#7b8aa5">本文实测</td><td>相对误差</td></tr>' +
      rows.map((r) => '<tr><td>' + r[0] + '</td><td style="text-align:left">' + r[1] +
        '</td><td style="text-align:left">' + r[2] + '</td><td>' + r[3] + '</td></tr>').join('') + '</table>';
  }

  /* ============================ 交互 ============================ */
  function bindRange(id, key, obj, render, digits) {
    const el = $(id), out = $(id + '-val');
    const sync = () => {
      obj[key] = parseFloat(el.value);
      out.textContent = obj[key].toFixed(digits === undefined ? 3 : digits);
      render();
    };
    el.addEventListener('input', sync);
    out.textContent = parseFloat(el.value).toFixed(digits === undefined ? 3 : digits);
  }

  function bindD2Drag() {
    const cv = plots.p2.cv;
    const local = (ev) => { const r = cv.getBoundingClientRect(); return [ev.clientX - r.left, ev.clientY - r.top]; };
    const hit = (px, py) => {
      for (let i = 0; i < D2.pts.length; i++) {
        const dx = plots.p2.X(D2.pts[i].re) - px, dy = plots.p2.Y(D2.pts[i].im) - py;
        if (dx * dx + dy * dy < 220) return i;
      }
      return -1;
    };
    cv.addEventListener('pointerdown', (ev) => {
      const p = local(ev), i = hit(p[0], p[1]);
      if (i < 0) return;
      D2.drag = i; cv.classList.add('dragging');
      cv.setPointerCapture(ev.pointerId);
      renderD2();
    });
    cv.addEventListener('pointermove', (ev) => {
      const p = local(ev);
      if (D2.drag >= 0) {
        const w = plots.p2.world;
        D2.pts[D2.drag].re = Math.max(w.xmin + 0.2, Math.min(w.xmax - 0.2, plots.p2.invX(p[0])));
        D2.pts[D2.drag].im = Math.max(w.ymin + 0.2, Math.min(w.ymax - 0.2, plots.p2.invY(p[1])));
        renderD2();
      } else {
        cv.style.cursor = hit(p[0], p[1]) >= 0 ? 'grab' : 'default';
      }
    });
    const up = () => { D2.drag = -1; cv.classList.remove('dragging'); renderD2(); };
    cv.addEventListener('pointerup', up);
    cv.addEventListener('pointercancel', up);
  }

  function boot() {
    ['p1a', 'p1b', 'p2', 'p3a', 'p3b', 'p4a', 'p4b'].forEach((id) => { plots[id] = new Plot(document.getElementById(id)); });

    bindRange('#eta1', 'eta', D1, renderD1);
    $('#rescale1').addEventListener('change', (ev) => { D1.rescale = ev.target.checked; renderD1(); });

    bindRange('#eta2', 'eta', D2, renderD2);
    $$('[data-preset]').forEach((btn) => btn.addEventListener('click', () => {
      D2.pts = D2_PRESETS[btn.dataset.preset].map((p) => ({ re: p.re, im: p.im }));
      $$('[data-preset]').forEach((b) => b.classList.toggle('active', b === btn));
      renderD2();
    }));
    bindD2Drag();

    bindRange('#eta3', 'eta', D3, renderD3, 2);
    bindRange('#beta3', 'beta', D3, renderD3, 2);
    $$('[data-d3]').forEach((btn) => btn.addEventListener('click', () => {
      const s = JSON.parse(btn.dataset.d3);
      $('#eta3').value = s.eta; $('#beta3').value = s.beta;
      D3.eta = s.eta; D3.beta = s.beta;
      $('#eta3-val').textContent = s.eta.toFixed(2);
      $('#beta3-val').textContent = s.beta.toFixed(2);
      renderD3();
    }));

    bindRange('#eta4', 'eta', D4, renderD4, 3);

    let raf = null;
    const renderAll = () => { renderD1(); renderD2(); renderD3(); renderD4(); };
    window.addEventListener('resize', () => {
      if (raf) cancelAnimationFrame(raf);
      raf = requestAnimationFrame(renderAll);
    });

    renderAll();
    renderSelfTest();
  }

  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', boot);
  else boot();
})();
