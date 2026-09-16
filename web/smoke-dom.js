/*!
 * smoke-dom.js — 无浏览器冒烟测试
 * ---------------------------------------------------------------------------
 * 用一个极小的 DOM 桩把 index.html 的元素结构喂给 stability-core.js + app.js，
 * 验证四个演示的渲染路径都能跑完、读数都写进了 DOM，并且交互（滑块 / 预设按钮 /
 * resize）重复渲染不抛异常。
 *
 * 运行： node smoke-dom.js
 *
 * 注意：这只是"渲染路径不炸"的验证，不能替代真机看效果 —— 布局与视觉仍需浏览器确认。
 */
'use strict';

/* ---------------- 最小 DOM 桩 ---------------- */
const noop = () => {};
function makeCtx(cv) {
  const ctx = {
    canvas: cv, fillStyle: '', strokeStyle: '', lineWidth: 1, font: '', textAlign: '', textBaseline: '',
    save: noop, restore: noop, beginPath: noop, moveTo: noop, lineTo: noop, stroke: noop, fill: noop,
    arc: noop, fillRect: noop, clearRect: noop, clip: noop, rect: noop,
    setTransform: noop, rotate: noop, translate: noop, fillText: noop, setLineDash: noop,
    measureText: () => ({ width: 10 }),
  };
  return ctx;
}
function makeCanvas(id) {
  const listeners = {};
  return {
    id, clientWidth: 352, clientHeight: 300, width: 0, height: 0, _listeners: listeners,
    getContext() { if (!this._ctx) this._ctx = makeCtx(this); return this._ctx; },
    addEventListener(t, f) { (listeners[t] = listeners[t] || []).push(f); },
    getBoundingClientRect() { return { left: 0, top: 0, width: this.clientWidth, height: this.clientHeight }; },
    setPointerCapture: noop,
    classList: { add: noop, remove: noop, toggle: noop },
    style: {},
  };
}
function makeEl(id, value) {
  const listeners = {};
  return {
    id, innerHTML: '', textContent: '', value: value === undefined ? '' : value,
    checked: true, dataset: {}, _listeners: listeners, style: {},
    classList: { add: noop, remove: noop, toggle: noop },
    addEventListener(t, f) { (listeners[t] = listeners[t] || []).push(f); },
  };
}
function fire(el, type, ev) {
  (el._listeners[type] || []).forEach((f) => f(Object.assign({ target: el }, ev || {})));
}

const canvases = ['p1a', 'p1b', 'p2', 'p3a', 'p3b', 'p4a', 'p4b'];
const els = {};
canvases.forEach((c) => { els[c] = makeCanvas(c); });
const defaults = { eta1: '0.15', eta2: '0.30', eta3: '0.60', beta3: '0.50', eta4: '0.150', rescale1: 'on' };
Object.keys(defaults).forEach((k) => { els[k] = makeEl(k, defaults[k]); });
['eta1-val', 'eta2-val', 'eta3-val', 'beta3-val', 'eta4-val',
  'd1-info', 'd1-tip', 'd2-table', 'd2-tip', 'd3-info', 'd3-tip',
  'd4-info', 'd4-tip', 'selftest'].forEach((k) => { els[k] = makeEl(k); });

const presets = ['real', 'complex', 'ill'].map((p) => Object.assign(makeEl('btn'), { dataset: { preset: p } }));
const d3btns = [
  { dataset: { d3: '{"eta":0.6,"beta":0.5}' } },
  { dataset: { d3: '{"eta":2.25,"beta":0.5}' } },
  { dataset: { d3: '{"eta":3.0,"beta":0.9}' } },
].map((b) => Object.assign(makeEl('btn'), b));

global.window = { devicePixelRatio: 1, addEventListener: noop };
global.requestAnimationFrame = () => 1;
global.cancelAnimationFrame = noop;
global.document = {
  readyState: 'complete',
  getElementById: (id) => els[id] || null,
  querySelector: (sel) => (sel.charAt(0) === '#' ? els[sel.slice(1)] || null : null),
  querySelectorAll: (sel) => (sel === '[data-preset]' ? presets : sel === '[data-d3]' ? d3btns : []),
  addEventListener: noop,
};

/* ---------------- 加载被测代码 ---------------- */
require('./stability-core.js');
const t0 = Date.now();
require('./app.js');
const bootMs = Date.now() - t0;

/* ---------------- 断言 ---------------- */
let pass = 0, fail = 0;
function ok(name, cond, extra) {
  cond ? pass++ : fail++;
  console.log((cond ? 'PASS  ' : 'FAIL  ') + name + (extra ? '   ' + extra : ''));
}
const has = (id, needle) => els[id].innerHTML.indexOf(needle) >= 0;

ok('boot() 无异常完成', true, '耗时 ' + bootMs + ' ms');
ok('① 读数表已填充（含 η_max）', has('d1-info', 'η_max') && has('d1-info', '谱半径'));
ok('① 结论文案已填充', els['d1-tip'].innerHTML.length > 40);
ok('② λ 判据表包含所有点', (els['d2-table'].innerHTML.match(/λ\d/g) || []).length >= 6);
ok('② 结论文案含圆盘判据', has('d2-tip', '圆盘'));
ok('③ 读数表含特征值与谱半径', has('d3-info', '谱半径') && has('d3-info', '特征值'));
ok('④ 读数表含稳态方差与噪声地板', has('d4-info', '稳态方差') && has('d4-info', '噪声地板'));
ok('⑤ 页内自检生成 ≥ 8 行', (els['selftest'].innerHTML.match(/<tr>/g) || []).length >= 8);
ok('⑤ 自检覆盖 A/B/C 三组', has('selftest', 'A ·') && has('selftest', 'B ·') && has('selftest', 'C ·'));

/* 交互：拖动滑块 → 重渲染；重复调用不应抛异常 */
const before = els['d1-info'].innerHTML;
els.eta1.value = '0.45';
fire(els.eta1, 'input');
ok('交互：① 滑块改变后读数更新', els['d1-info'].innerHTML !== before && has('d1-info', '0.450'));
ok('交互：① η < η_max 判为收敛', has('d1-info', '离散收敛'));

els.eta1.value = '0.6';
fire(els.eta1, 'input');
ok('交互：① η > η_max 判为发散', has('d1-info', '离散发散'));

els.rescale1.checked = false;
fire(els.rescale1, 'change');
ok('交互：① 关闭重标定后文案切换', els['d1-tip'].innerHTML.indexOf('时间轴') >= 0);

fire(presets[1], 'click');
ok('交互：② 预设（复特征值）生效', (els['d2-table'].innerHTML.match(/λ\d/g) || []).length === 4);

fire(d3btns[1], 'click');
ok('交互：③ GAP 反例预设生效', has('d3-info', '2.25') && has('d3-tip', '缝隙'));

els.eta2.value = '0.6';
fire(els.eta2, 'input');
ok('交互：② η 调大后出现缝隙提示', els['d2-tip'].innerHTML.length > 40);

/* canvas 指针事件：拖拽 λ 点（② 演示） */
const pc = els.p2;
fire(pc, 'pointerdown', { clientX: 120, clientY: 150, pointerId: 1 });
fire(pc, 'pointermove', { clientX: 160, clientY: 160, pointerId: 1 });
fire(pc, 'pointerup', { clientX: 160, clientY: 160, pointerId: 1 });
ok('交互：② 拖拽事件路径无异常', true);

console.log('='.repeat(60));
console.log(pass + ' passed, ' + fail + ' failed');
process.exit(fail ? 1 : 0);
