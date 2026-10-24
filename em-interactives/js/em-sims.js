/*
 * Interactive electron-microscopy demos: FIB-SEM tomography, 4D-STEM and EELS
 * for tofunmi-ogunfunmi.github.io
 *
 * Written from scratch for this site. The visual approach (oblique 3D scene of
 * probe, sample and detector; animated "pillar" cards) is inspired by the
 * Colin Ophus Lab interactive 4D-STEM demo: https://colab.stanford.edu/interactive
 *
 * These are simplified teaching models, not quantitative simulations.
 */
(function () {
'use strict';

/* ------------------------------------------------------------------ utils */
const TAU = Math.PI * 2, D2R = Math.PI / 180;
const clamp = (v, a, b) => (v < a ? a : v > b ? b : v);
const REDUCED = !!(window.matchMedia && matchMedia('(prefers-reduced-motion: reduce)').matches);
const FONT = '"DM Sans", system-ui, sans-serif';
const COL = { panel: '#0c0f14', t1: '#f0ede8', t2: '#9a968f', acc: '#58b947', blu: '#4a9eff', amb: '#f5b83d', red: '#ff6b5a' };

function rng32(a) { return function () { a |= 0; a = a + 0x6D2B79F5 | 0; let t = Math.imul(a ^ a >>> 15, 1 | a); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; }; }
function gauss(r) { let u = 0; while (u === 0) u = r(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(TAU * r()); }
function erf(x) { const s = x < 0 ? -1 : 1; x = Math.abs(x); const t = 1 / (1 + 0.3275911 * x); const y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x); return s * y; }
function hsv(h, s, v) { h = ((h % 1) + 1) % 1; const i = Math.floor(h * 6), f = h * 6 - i, p = v * (1 - s), q = v * (1 - f * s), t = v * (1 - (1 - f) * s); switch (i % 6) { case 0: return [v, t, p]; case 1: return [q, v, p]; case 2: return [p, v, t]; case 3: return [p, q, v]; case 4: return [t, p, v]; default: return [v, p, q]; } }
function hash3(i, j, k) { let h = (Math.imul(i, 374761393) + Math.imul(j, 668265263) + Math.imul(k, 1274126177)) | 0; h = Math.imul(h ^ (h >>> 13), 1274126177); h ^= h >>> 16; return (h >>> 0) / 4294967296; }
function vnoise3(x, y, z) {
  const i = Math.floor(x), j = Math.floor(y), k = Math.floor(z);
  let fx = x - i, fy = y - j, fz = z - k; fx = fx * fx * (3 - 2 * fx); fy = fy * fy * (3 - 2 * fy); fz = fz * fz * (3 - 2 * fz);
  const L = (a, b, t) => a + (b - a) * t;
  const c00 = L(hash3(i, j, k), hash3(i + 1, j, k), fx), c10 = L(hash3(i, j + 1, k), hash3(i + 1, j + 1, k), fx);
  const c01 = L(hash3(i, j, k + 1), hash3(i + 1, j, k + 1), fx), c11 = L(hash3(i, j + 1, k + 1), hash3(i + 1, j + 1, k + 1), fx);
  return L(L(c00, c10, fy), L(c01, c11, fy), fz);
}
const SUP = { '0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹', '-': '⁻' };
function sci(x, d) { if (!isFinite(x)) return '–'; if (Math.abs(x) < 1e4) return x.toFixed(d == null ? 0 : d); const e = Math.floor(Math.log10(Math.abs(x))); return (x / Math.pow(10, e)).toFixed(1) + '×10' + String(e).split('').map(c => SUP[c]).join(''); }

/* inferno-like colormap (matplotlib "inferno" control points) */
const LUT = (function () {
  const s = [[0, 0, 4], [31, 12, 72], [85, 15, 109], [136, 34, 106], [186, 54, 85], [227, 89, 51], [249, 140, 10], [249, 201, 50], [252, 255, 164]];
  const L = new Uint8ClampedArray(768);
  for (let i = 0; i < 256; i++) { const t = i / 255 * (s.length - 1), j = Math.min(s.length - 2, Math.floor(t)), f = t - j; for (let c = 0; c < 3; c++) L[i * 3 + c] = s[j][c] + (s[j + 1][c] - s[j][c]) * f; }
  return L;
})();

/* canvas helpers: fixed logical size, sharp on hi-dpi, scales with CSS width */
function makeCanvas(cv, w, h) {
  const dpr = Math.min(2, window.devicePixelRatio || 1);
  cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr);
  const ctx = cv.getContext('2d');
  return { cv, ctx, w, h, dpr, reset() { ctx.setTransform(dpr, 0, 0, dpr, 0, 0); } };
}
function offscreen(w, h) { const c = document.createElement('canvas'); c.width = w; c.height = h; return c; }
function pos(V, e) { const r = V.cv.getBoundingClientRect(); return [(e.clientX - r.left) / r.width * V.w, (e.clientY - r.top) / r.height * V.h]; }
/* draw an image onto the parallelogram p00 (top-left) -> p10 (top-right), p01 (bottom-left) */
function quad(V, img, p00, p10, p01, alpha) {
  const { ctx, dpr } = V, w = img.width, h = img.height;
  ctx.save(); if (alpha != null) ctx.globalAlpha = alpha;
  ctx.setTransform(dpr * (p10[0] - p00[0]) / w, dpr * (p10[1] - p00[1]) / w, dpr * (p01[0] - p00[0]) / h, dpr * (p01[1] - p00[1]) / h, dpr * p00[0], dpr * p00[1]);
  ctx.drawImage(img, 0, 0); ctx.restore();
}
function poly(ctx, p) { ctx.beginPath(); ctx.moveTo(p[0][0], p[0][1]); for (let i = 1; i < p.length; i++) ctx.lineTo(p[i][0], p[i][1]); ctx.closePath(); }
function hull(pts) {
  const p = pts.slice().sort((a, b) => a[0] - b[0] || a[1] - b[1]);
  const cr = (o, a, b) => (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
  const lo = [], up = [];
  for (const q of p) { while (lo.length >= 2 && cr(lo[lo.length - 2], lo[lo.length - 1], q) <= 0) lo.pop(); lo.push(q); }
  for (let i = p.length - 1; i >= 0; i--) { const q = p[i]; while (up.length >= 2 && cr(up[up.length - 2], up[up.length - 1], q) <= 0) up.pop(); up.push(q); }
  up.pop(); lo.pop(); return lo.concat(up);
}
function txt(ctx, s, x, y, o) { o = o || {}; ctx.font = (o.w || 500) + ' ' + (o.size || 11) + 'px ' + FONT; ctx.fillStyle = o.color || COL.t2; ctx.textAlign = o.align || 'left'; ctx.textBaseline = o.base || 'middle'; ctx.fillText(s, x, y); }
function scaleBar(ctx, x1, x2, y, label) { ctx.strokeStyle = '#fff'; ctx.lineWidth = 2; ctx.beginPath(); ctx.moveTo(x1, y); ctx.lineTo(x2, y); ctx.stroke(); txt(ctx, label, (x1 + x2) / 2, y - 9, { align: 'center', color: '#fff', size: 10, w: 600 }); }
function circ3(cx, cy, z, r, n) { const o = []; for (let i = 0; i < n; i++) { const a = i / n * TAU; o.push([cx + r * Math.cos(a), cy + r * Math.sin(a), z]); } return o; }

/* oblique parallel-projection camera; x right, y away from viewer, z up */
function Cam(yaw, elev) { this.s = 1; this.cx = 0; this.cy = 0; this.set(yaw, elev); }
Cam.prototype.set = function (yaw, elev) { this.yaw = yaw; this.elev = elev; this.cY = Math.cos(yaw * D2R); this.sY = Math.sin(yaw * D2R); this.cE = Math.cos(elev * D2R); this.sE = Math.sin(elev * D2R); };
Cam.prototype.p = function (x, y, z) { const xr = x * this.cY - y * this.sY, yr = x * this.sY + y * this.cY; return [this.cx + this.s * xr, this.cy - this.s * (z * this.cE + yr * this.sE)]; };
Cam.prototype.depth = function (x, y, z) { return (x * this.sY + y * this.cY) * this.cE - z * this.sE; };
Cam.prototype.unproject = function (X, Y, z) { const xr = (X - this.cx) / this.s, yr = ((this.cy - Y) / this.s - z * this.cE) / this.sE; return [xr * this.cY + yr * this.sY, -xr * this.sY + yr * this.cY]; };
Cam.prototype.fit = function (pts, W, H, pad) {
  this.s = 1; this.cx = 0; this.cy = 0; let a = 1e9, b = -1e9, c = 1e9, d = -1e9;
  for (const q of pts) { const r = this.p(q[0], q[1], q[2]); a = Math.min(a, r[0]); b = Math.max(b, r[0]); c = Math.min(c, r[1]); d = Math.max(d, r[1]); }
  const s = Math.min((W - 2 * pad) / (b - a), (H - 2 * pad) / (d - c)); this.s = s; this.cx = W / 2 - s * (a + b) / 2; this.cy = H / 2 - s * (c + d) / 2;
};
function cone(ctx, cam, a, b, stroke) {
  const pts = [];
  for (const n of [a, b]) { if (n[3] < 0.05) pts.push(cam.p(n[0], n[1], n[2])); else for (const q of circ3(n[0], n[1], n[2], n[3], 24)) pts.push(cam.p(q[0], q[1], q[2])); }
  poly(ctx, hull(pts)); ctx.fill(); if (stroke) ctx.stroke();
}

/* radix-2 complex FFT */
function FFT(n) {
  this.n = n; const lv = Math.round(Math.log2(n)); this.rev = new Uint32Array(n);
  for (let i = 0; i < n; i++) { let r = 0; for (let b = 0; b < lv; b++) r |= ((i >> b) & 1) << (lv - 1 - b); this.rev[i] = r; }
  this.cos = new Float64Array(n / 2); this.sin = new Float64Array(n / 2);
  for (let i = 0; i < n / 2; i++) { this.cos[i] = Math.cos(TAU * i / n); this.sin[i] = Math.sin(TAU * i / n); }
  this.tr = new Float64Array(n); this.ti = new Float64Array(n);
}
FFT.prototype.f1 = function (re, im, inv) {
  const n = this.n, rev = this.rev, C = this.cos, S = this.sin;
  for (let i = 0; i < n; i++) { const j = rev[i]; if (j > i) { let t = re[i]; re[i] = re[j]; re[j] = t; t = im[i]; im[i] = im[j]; im[j] = t; } }
  for (let size = 2; size <= n; size <<= 1) {
    const half = size >> 1, step = n / size;
    for (let i = 0; i < n; i += size) for (let j = 0, k = 0; j < half; j++, k += step) {
      const wr = C[k], wi = inv ? S[k] : -S[k], a = i + j, b = a + half;
      const xr = re[b] * wr - im[b] * wi, xi = re[b] * wi + im[b] * wr;
      re[b] = re[a] - xr; im[b] = im[a] - xi; re[a] += xr; im[a] += xi;
    }
  }
};
FFT.prototype.f2 = function (re, im, inv) {
  const n = this.n, tr = this.tr, ti = this.ti;
  for (let y = 0; y < n; y++) { const o = y * n; for (let x = 0; x < n; x++) { tr[x] = re[o + x]; ti[x] = im[o + x]; } this.f1(tr, ti, inv); for (let x = 0; x < n; x++) { re[o + x] = tr[x]; im[o + x] = ti[x]; } }
  for (let x = 0; x < n; x++) { for (let y = 0; y < n; y++) { tr[y] = re[y * n + x]; ti[y] = im[y * n + x]; } this.f1(tr, ti, inv); for (let y = 0; y < n; y++) { re[y * n + x] = tr[y]; im[y * n + x] = ti[y]; } }
};

/* only animate while on screen (hidden pages / scrolled away cost nothing) */
function whenVisible(el, frame) {
  let on = false, id = 0, last = 0;
  const loop = t => { if (!on) return; const dt = last ? Math.min(0.1, (t - last) / 1000) : 0.016; last = t; try { frame(dt); } catch (e) { on = false; console.error(e); return; } id = requestAnimationFrame(loop); };
  if (!('IntersectionObserver' in window)) { on = true; id = requestAnimationFrame(loop); return; }
  new IntersectionObserver(es => { for (const e of es) { if (e.isIntersecting && !on) { on = true; last = 0; id = requestAnimationFrame(loop); } else if (!e.isIntersecting && on) { on = false; cancelAnimationFrame(id); } } }, { rootMargin: '80px' }).observe(el);
}
function drag(V, fn, opt) {
  let down = false; const c = V.cv;
  c.addEventListener('pointerdown', e => { down = true; try { c.setPointerCapture(e.pointerId); } catch (_) {} fn(pos(V, e), 'down', e); });
  c.addEventListener('pointermove', e => { if (down) fn(pos(V, e), 'move', e); });
  const up = e => { if (down) { down = false; fn(pos(V, e), 'up', e); } };
  c.addEventListener('pointerup', up); c.addEventListener('pointercancel', up);
}

/* tiny DOM builders for the controls */
function h(tag, attrs, kids) {
  const e = document.createElement(tag);
  if (attrs) for (const k in attrs) { const v = attrs[k]; if (v == null) continue; if (k === 'class') e.className = v; else if (k === 'text') e.textContent = v; else if (k.slice(0, 2) === 'on') e.addEventListener(k.slice(2), v); else e.setAttribute(k, v); }
  if (kids != null) (Array.isArray(kids) ? kids : [kids]).forEach(c => { if (c != null) e.append(c); });
  return e;
}
let UID = 0;
function slider(o) {
  const id = 'emsim-' + (++UID), out = h('output', { for: id });
  const inp = h('input', { type: 'range', id, min: o.min, max: o.max, step: o.step, value: o.value });
  const upd = () => { out.textContent = o.fmt(+inp.value); };
  inp.addEventListener('input', () => { upd(); o.on(+inp.value); }); upd();
  const w = h('div', { class: 'ctrl' }, [h('label', { for: id }, [h('span', { text: o.label }), out]), inp]);
  w.set = v => { inp.value = v; upd(); }; return w;
}
function seg(o) {
  const w = h('div', { class: 'ctrl' }), g = h('div', { class: 'seg', role: 'group', 'aria-label': o.label }), btns = {};
  w.append(h('span', { class: 'ctrl-l', text: o.label }));
  o.options.forEach(([v, t]) => { const b = h('button', { type: 'button', 'aria-pressed': String(v === o.value), text: t, onclick: () => { w.set(v); o.on(v); } }); btns[v] = b; g.append(b); });
  w.append(g); w.set = v => { for (const k in btns) btns[k].setAttribute('aria-pressed', String(k === String(v))); }; return w;
}
function btn(text, on, cls) { return h('button', { type: 'button', class: 'sbtn' + (cls ? ' ' + cls : ''), text, onclick: on }); }
function fail(root, e) { console.error(e); root.append(h('p', { class: 'hint', text: 'This demo could not start in your browser.' })); }

/* ================================================================ 4D-STEM */
const ST = { px: 0.25, N: 128, FW: 384, LAMBDA: 0.0197 }; // 300 kV electrons, 0.25 Å pixels
ST.L = ST.FW * ST.px; ST.dk = 1 / (ST.N * ST.px); ST.H = ST.L / 2; ST.S0 = ST.N / 2; ST.S1 = ST.FW - ST.N / 2;
const GRAIN_RGB = [[255, 122, 89], [63, 208, 224], [185, 168, 216], [245, 184, 61], [120, 214, 110]];
let STEM_SAMPLE = null;
function stemSample() {
  if (STEM_SAMPLE) return STEM_SAMPLE;
  const { px, FW, L } = ST, r = rng32(11);
  // four crystalline grains (three lattices) and one amorphous pocket
  const seeds = [{ x: 20, y: 20, t: 'sq', a: 2.0, ang: 15 }, { x: 77, y: 17, t: 'hex', a: 2.45, ang: 22.5 }, { x: 52, y: 54, t: 'am' }, { x: 17, y: 79, t: 'rect', a: 2.1, b: 2.9, ang: -22.5 }, { x: 81, y: 80, t: 'sq', a: 2.0, ang: 30 }];
  const gAt = (x, y) => { const w = 2.2 * Math.sin(x * 0.13 + y * 0.07) + 1.8 * Math.sin(y * 0.19 - x * 0.06 + 1.3); let b = 0, bd = 1e9; for (let i = 0; i < seeds.length; i++) { const s = seeds[i]; const d = Math.hypot(x - s.x, y - s.y) + (i % 2 ? w : -w) + (s.t === 'am' ? 5 : 0); if (d < bd) { bd = d; b = i; } } return b; };
  const atoms = [];
  seeds.forEach((s, g) => {
    if (s.t === 'am') { // random close packing by dart throwing
      const cell = 2, md2 = 1.9 * 1.9, grid = new Map();
      for (let k = 0; k < 40000; k++) {
        const x = r() * L, y = r() * L; if (gAt(x, y) !== g) continue;
        const ci = Math.floor(x / cell), cj = Math.floor(y / cell); let ok = true;
        for (let di = -1; di <= 1 && ok; di++) for (let dj = -1; dj <= 1 && ok; dj++) { const arr = grid.get((ci + di) * 4096 + cj + dj); if (arr) for (let q = 0; q < arr.length; q += 2) { const dx = arr[q] - x, dy = arr[q + 1] - y; if (dx * dx + dy * dy < md2) { ok = false; break; } } }
        if (!ok) continue; const key = ci * 4096 + cj; let arr = grid.get(key); if (!arr) { arr = []; grid.set(key, arr); } arr.push(x, y); atoms.push(x, y, 1.05, g);
      }
    } else {
      const an = s.ang * D2R, c = Math.cos(an), sn = Math.sin(an), a = s.a, b = s.b || s.a;
      const v1 = [a, 0], v2 = s.t === 'hex' ? [a / 2, a * Math.sqrt(3) / 2] : [0, b];
      const r1 = [v1[0] * c - v1[1] * sn, v1[0] * sn + v1[1] * c], r2 = [v2[0] * c - v2[1] * sn, v2[0] * sn + v2[1] * c];
      const basis = s.t === 'rect' ? [[0, 0, 1.0], [0.5, 0.5, 0.55]] : [[0, 0, s.t === 'hex' ? 1.1 : 1.0]];
      const M = Math.ceil(L * 1.5 / Math.min(a, b));
      for (let i = -M; i <= M; i++) for (let j = -M; j <= M; j++) for (const bs of basis) {
        const fi = i + bs[0], fj = j + bs[1], x = s.x + fi * r1[0] + fj * r2[0], y = s.y + fi * r1[1] + fj * r2[1];
        if (x < 0 || y < 0 || x >= L || y >= L || gAt(x, y) !== g) continue; atoms.push(x, y, bs[2], g);
      }
    }
  });
  const V = new Float32Array(FW * FW), sig = 0.32 / px, R = Math.ceil(3 * sig), i2s = 1 / (2 * sig * sig);
  for (let k = 0; k < atoms.length; k += 4) {
    const ax = atoms[k] / px, ay = atoms[k + 1] / px, w = atoms[k + 2];
    const x0 = Math.max(0, Math.ceil(ax - R)), x1 = Math.min(FW - 1, Math.floor(ax + R)), y0 = Math.max(0, Math.ceil(ay - R)), y1 = Math.min(FW - 1, Math.floor(ay + R));
    for (let y = y0; y <= y1; y++) { const dy = y - ay; for (let x = x0; x <= x1; x++) { const dx = x - ax; V[y * FW + x] += w * Math.exp(-(dx * dx + dy * dy) * i2s); } }
  }
  let vmax = 0; for (let i = 0; i < V.length; i++) if (V[i] > vmax) vmax = V[i];
  const tRe = new Float32Array(FW * FW), tIm = new Float32Array(FW * FW), tex = offscreen(FW, FW), tctx = tex.getContext('2d'), img = tctx.createImageData(FW, FW), d = img.data;
  for (let y = 0; y < FW; y++) for (let x = 0; x < FW; x++) {
    const i = y * FW + x, v = V[i] / vmax, ph = 1.1 * v; tRe[i] = Math.cos(ph); tIm[i] = Math.sin(ph);
    const c = GRAIN_RGB[gAt(x * px, y * px)], q = Math.pow(clamp(v * 1.3, 0, 1), 0.85);
    d[i * 4] = 14 + (c[0] - 14) * q; d[i * 4 + 1] = 17 + (c[1] - 17) * q; d[i * 4 + 2] = 24 + (c[2] - 24) * q; d[i * 4 + 3] = 255;
  }
  tctx.putImageData(img, 0, 0);
  return (STEM_SAMPLE = { tRe, tIm, tex });
}

function StemEngine(alpha) {
  const N = ST.N; this.S = stemSample(); this.fft = new FFT(N);
  this.pRe = new Float64Array(N * N); this.pIm = new Float64Array(N * N); this.wRe = new Float64Array(N * N); this.wIm = new Float64Array(N * N);
  this.I = new Float32Array(N * N); this.alpha = alpha; this.df = 0; this.makeProbe();
}
/* probe = inverse FFT of a (soft-edged) aperture with a defocus phase */
StemEngine.prototype.makeProbe = function () {
  const N = ST.N, dk = ST.dk, h = N / 2, re = this.wRe, im = this.wIm, a = this.alpha, df = this.df;
  for (let y = 0; y < N; y++) { const ky = (y < h ? y : y - N) * dk; for (let x = 0; x < N; x++) {
    const kx = (x < h ? x : x - N) * dk, k2 = kx * kx + ky * ky, A = clamp((a - Math.sqrt(k2)) / dk + 0.5, 0, 1), chi = Math.PI * ST.LAMBDA * df * k2;
    re[y * N + x] = A * Math.cos(chi); im[y * N + x] = -A * Math.sin(chi);
  } }
  this.fft.f2(re, im, true);
  let s = 0; for (let i = 0; i < N * N; i++) s += re[i] * re[i] + im[i] * im[i];
  const nr = 1 / Math.sqrt(s);
  for (let y = 0; y < N; y++) for (let x = 0; x < N; x++) { const src = ((y + h) & (N - 1)) * N + ((x + h) & (N - 1)); this.pRe[y * N + x] = re[src] * nr; this.pIm[y * N + x] = im[src] * nr; }
  // annulus used to estimate in-plane orientation from Bragg spots (12-fold harmonic covers 4- and 6-fold patterns)
  const r0 = (a + 0.07) / dk, r1 = Math.min(1.25 / dk, h - 2), idx = [], cs = [], sn = [];
  for (let y = 0; y < N; y++) for (let x = 0; x < N; x++) { const dx = x - h, dy = y - h, rr = Math.sqrt(dx * dx + dy * dy); if (rr < r0 || rr > r1) continue; const g = 12 * Math.atan2(dy, dx); idx.push(y * N + x); cs.push(Math.cos(g)); sn.push(Math.sin(g)); }
  this.oIdx = Int32Array.from(idx); this.oCos = Float32Array.from(cs); this.oSin = Float32Array.from(sn); this.oR = [r0 * dk, r1 * dk];
};
/* exit wave = probe x transmission function; diffraction pattern = |FFT|^2 (fft-shifted) */
StemEngine.prototype.diffract = function (cx, cy) {
  const N = ST.N, FW = ST.FW, h = N / 2, re = this.wRe, im = this.wIm, pR = this.pRe, pI = this.pIm, tR = this.S.tRe, tI = this.S.tIm, out = this.I;
  const x0 = cx - h, y0 = cy - h;
  for (let y = 0; y < N; y++) { const fo = (y0 + y) * FW + x0, o = y * N; for (let x = 0; x < N; x++) { const a = pR[o + x], b = pI[o + x], c = tR[fo + x], d = tI[fo + x]; re[o + x] = a * c - b * d; im[o + x] = a * d + b * c; } }
  this.fft.f2(re, im, false);
  for (let y = 0; y < N; y++) { const sy = ((y + h) & (N - 1)) * N; for (let x = 0; x < N; x++) { const i = y * N + x; out[sy + ((x + h) & (N - 1))] = re[i] * re[i] + im[i] * im[i]; } }
  return out;
};
StemEngine.prototype.orient = function (I) {
  const idx = this.oIdx, C = this.oCos, S = this.oSin; let zr = 0, zi = 0, s = 0;
  for (let j = 0; j < idx.length; j++) { const v = I[idx[j]], w = v * v; zr += w * C[j]; zi += w * S[j]; s += w; }
  return [(Math.atan2(zi, zr) / TAU + 1) % 1, s > 0 ? Math.hypot(zr, zi) / s : 0];
};
function scanPos(i, SG) { const sx = i % SG, sy = (i / SG) | 0, st = (ST.S1 - ST.S0) / (SG - 1); return [Math.round(ST.S0 + sx * st), Math.round(ST.S0 + sy * st)]; }
function paintDP(I, img, gamma) {
  let mx = 0; for (let i = 0; i < I.length; i++) if (I[i] > mx) mx = I[i];
  const d = img.data, inv = mx > 0 ? 1 / mx : 0;
  for (let i = 0; i < I.length; i++) { const li = Math.min(255, (Math.pow(I[i] * inv, gamma) * 255) | 0) * 3; d[i * 4] = LUT[li]; d[i * 4 + 1] = LUT[li + 1]; d[i * 4 + 2] = LUT[li + 2]; d[i * 4 + 3] = 255; }
}
const ORI = c => clamp((c - 0.2) / 0.25, 0, 1);
function oriRGB(hue, coh) { return hsv(hue, 0.78, 0.12 + 0.88 * ORI(coh)); }
function composeStem(cctx, tex, vCv, SG) {
  const st = (ST.S1 - ST.S0) / (SG - 1);
  cctx.clearRect(0, 0, ST.FW, ST.FW); cctx.drawImage(tex, 0, 0);
  cctx.imageSmoothingEnabled = false; cctx.globalAlpha = 0.93;
  cctx.drawImage(vCv, ST.S0 - st / 2, ST.S0 - st / 2, ST.S1 - ST.S0 + st, ST.S1 - ST.S0 + st);
  cctx.globalAlpha = 1; cctx.imageSmoothingEnabled = true;
}
const STEM_GEO = { zIn: 60, zDet: -88, Dh: 46, zB: -8 }, STEM_GEO_SMALL = { zIn: 44, zDet: -62, Dh: 40, zB: -6 };
function stemCamFit(cam, w, h, pad, G) {
  const H = ST.H, pts = []; G = G || STEM_GEO;
  for (const x of [-H, H]) for (const y of [-H, H]) pts.push([x, y, 0], [x, y, G.zB]);
  for (const x of [-G.Dh, G.Dh]) for (const y of [-G.Dh, G.Dh]) pts.push([x, y, G.zDet]);
  for (const q of circ3(0, 0, G.zIn, 30, 16)) pts.push(q);
  pts.push([0, 0, G.zIn + 34]); cam.fit(pts, w, h, pad);
}
function drawStemScene(V, eng, comp, dpCv, probe, labels) {
  const { ctx, w, h } = V, cam = V.cam, H = ST.H, G = V.geo || STEM_GEO, P = (x, y, z) => cam.p(x, y, z);
  V.reset(); ctx.fillStyle = COL.panel; ctx.fillRect(0, 0, w, h);
  const M = V.map, wp = M ? M.world(probe) : [probe[0] * ST.px - H, H - probe[1] * ST.px], wx = wp[0], wy = wp[1], kmax = M ? M.kmax : 2;
  const Rin = 5 + eng.alpha * 26, zc = clamp(eng.df * 0.05, -30, 30), k = Rin / (G.zIn - zc);
  // detector with the live diffraction pattern
  const d00 = P(-G.Dh, G.Dh, G.zDet), d10 = P(G.Dh, G.Dh, G.zDet), d01 = P(-G.Dh, -G.Dh, G.zDet), d11 = P(G.Dh, -G.Dh, G.zDet);
  ctx.fillStyle = '#000'; poly(ctx, [d00, d10, d11, d01]); ctx.fill();
  quad(V, dpCv, d00, d10, d01);
  ctx.strokeStyle = 'rgba(255,255,255,.18)'; ctx.lineWidth = 1; poly(ctx, [d00, d10, d11, d01]); ctx.stroke();
  // transmitted beam below the sample
  const lo = [[wx, wy, G.zB, Math.max(0.5, Math.abs(G.zB - zc) * k)]];
  if (zc < G.zB - 0.5) { const f = (zc - G.zB) / (G.zDet - G.zB); lo.push([wx * (1 - f), wy * (1 - f), zc, 0]); }
  lo.push([0, 0, G.zDet, Math.max(1, eng.alpha / kmax * G.Dh)]);
  ctx.fillStyle = 'rgba(88,185,71,.12)'; for (let i = 0; i < lo.length - 1; i++) cone(ctx, cam, lo[i], lo[i + 1]);
  // sample slab
  const T = [[-H, H], [H, H], [H, -H], [-H, -H]], sides = [];
  for (let i = 0; i < 4; i++) { const a = T[i], b = T[(i + 1) % 4]; sides.push({ q: [P(a[0], a[1], 0), P(b[0], b[1], 0), P(b[0], b[1], G.zB), P(a[0], a[1], G.zB)], d: cam.depth((a[0] + b[0]) / 2, (a[1] + b[1]) / 2, G.zB / 2), i }); }
  sides.sort((u, v) => v.d - u.d);
  for (const s of sides) { ctx.fillStyle = s.i % 2 ? '#121722' : '#1a2130'; poly(ctx, s.q); ctx.fill(); ctx.strokeStyle = 'rgba(255,255,255,.07)'; ctx.stroke(); }
  const t00 = P(-H, H, 0), t10 = P(H, H, 0), t01 = P(-H, -H, 0), t11 = P(H, -H, 0);
  quad(V, comp, t00, t10, t01); ctx.strokeStyle = 'rgba(255,255,255,.2)'; poly(ctx, [t00, t10, t11, t01]); ctx.stroke();
  const sr = M ? M.scanH : ST.S0 * ST.px - H; // scan field outline
  if (sr) { ctx.setLineDash([3, 4]); ctx.strokeStyle = 'rgba(255,255,255,.28)'; poly(ctx, [P(sr, -sr, 0), P(-sr, -sr, 0), P(-sr, sr, 0), P(sr, sr, 0)]); ctx.stroke(); ctx.setLineDash([]); }
  // converging probe above the sample
  const up = [[wx, wy, G.zIn, Rin]]; if (zc > 0.5) up.push([wx, wy, zc, 0]); up.push([wx, wy, 0, Math.max(0.5, Math.abs(zc) * k)]);
  ctx.fillStyle = 'rgba(88,185,71,.2)'; ctx.strokeStyle = 'rgba(130,225,110,.5)'; ctx.lineWidth = 1;
  for (let i = 0; i < up.length - 1; i++) cone(ctx, cam, up[i], up[i + 1], true);
  ctx.strokeStyle = 'rgba(170,245,150,.75)'; poly(ctx, circ3(wx, wy, G.zIn, Rin, 40).map(q => P(q[0], q[1], q[2]))); ctx.stroke();
  const sp = P(wx, wy, 0); ctx.save(); ctx.shadowColor = '#7dff6a'; ctx.shadowBlur = 14; ctx.fillStyle = '#eaffe4'; ctx.beginPath(); ctx.arc(sp[0], sp[1], 2.6, 0, TAU); ctx.fill(); ctx.restore();
  if (labels) {
    const top = P(wx, wy, G.zIn); txt(ctx, 'electron probe', top[0], top[1] - Rin * cam.s * cam.sE - 10, { align: 'center', color: COL.t1 });
    const rc = [t00, t10, t01, t11].reduce((a, b) => (b[0] > a[0] ? b : a)); txt(ctx, 'sample', rc[0] + 8, rc[1] + 4);
    const dc = [d00, d10, d01, d11].reduce((a, b) => (b[0] > a[0] ? b : a)); txt(ctx, 'detector · diffraction pattern', dc[0] - 6, dc[1] + 26, { align: 'right' });
  }
}

/* pillar card: probe rasters across the sample building a crystal-orientation map */
function stemPillar(cv) {
  const V = makeCanvas(cv, 360, 360); V.cam = new Cam(-35, 32); V.geo = STEM_GEO_SMALL; let s = null; const SG = 32;
  function init() {
    const eng = new StemEngine(0.13), vCv = offscreen(SG, SG), vctx = vCv.getContext('2d'), comp = offscreen(ST.FW, ST.FW), dpCv = offscreen(ST.N, ST.N), dctx = dpCv.getContext('2d');
    stemCamFit(V.cam, 360, 360, 10, V.geo);
    s = { eng, vCv, vctx, vImg: vctx.createImageData(SG, SG), comp, cctx: comp.getContext('2d'), dpCv, dctx, dImg: dctx.createImageData(ST.N, ST.N), i: 0, hold: 0, probe: [ST.S0, ST.S0] };
  }
  function step(n) {
    const d = s.vImg.data;
    for (let k = 0; k < n && s.i < SG * SG; k++) { const p = scanPos(s.i, SG); s.eng.diffract(p[0], p[1]); const o = s.eng.orient(s.eng.I), c = oriRGB(o[0], o[1]); d[s.i * 4] = c[0] * 255; d[s.i * 4 + 1] = c[1] * 255; d[s.i * 4 + 2] = c[2] * 255; d[s.i * 4 + 3] = 245; s.probe = p; s.i++; }
  }
  function render() { s.vctx.putImageData(s.vImg, 0, 0); composeStem(s.cctx, s.eng.S.tex, s.vCv, SG); paintDP(s.eng.I, s.dImg, 0.28); s.dctx.putImageData(s.dImg, 0, 0); drawStemScene(V, s.eng, s.comp, s.dpCv, s.probe, false); }
  whenVisible(cv, () => {
    if (!s) init();
    if (REDUCED) { if (!s.done) { step(SG * 20); render(); s.done = true; } return; }
    if (s.hold > 0) { if (--s.hold === 0) { s.i = 0; s.vImg.data.fill(0); } } else { step(3); if (s.i >= SG * SG) s.hold = 150; }
    render();
  });
}

function stemWidget(root) {
  const SG = 40, NB = 40, NN = NB * NB;
  const sceneCv = h('canvas', { class: 'scene', role: 'img', 'aria-label': '3D view of a 4D-STEM experiment: a converging electron probe above a polycrystalline sample and the diffraction pattern on the detector below. Drag to move the probe.' });
  const aCv = h('canvas', { tabindex: '0', 'aria-label': 'Virtual image of the sample. Drag, or use the arrow keys, to move the probe.' });
  const bCv = h('canvas', { 'aria-label': 'Diffraction pattern at the probe position. Click a Bragg spot to make a dark-field image from it.' });
  const SC = makeCanvas(sceneCv, 760, 420), PA = makeCanvas(aCv, 360, 360), PB = makeCanvas(bCv, 360, 360);
  SC.cam = new Cam(-32, 28); stemCamFit(SC.cam, 760, 420, 16);
  const status = h('span', { class: 'hint', 'aria-live': 'polite' }), aTitle = h('span');
  let eng = null, cube, oriH, oriC, vals, acquired = 0, scanning = false, dirty = true, mode = 'adf', view = 'vimg', gamma = 0.3, probe = [96, 96], home = [96, 96], vdf = null, maskIdx = new Int32Array(0);
  const vdfR = 0.085, vCv = offscreen(SG, SG), vctx = vCv.getContext('2d'), vImg = vctx.createImageData(SG, SG);
  const dpCv = offscreen(ST.N, ST.N), dctx = dpCv.getContext('2d'), dImg = dctx.createImageData(ST.N, ST.N), comp = offscreen(ST.FW, ST.FW), cctx = comp.getContext('2d');
  const stp = (ST.S1 - ST.S0) / (SG - 1), R0 = ST.S0 - stp / 2, RS = ST.S1 - ST.S0 + stp;
  const adfRange = () => { const o = 1.22; return [Math.min(Math.max(eng.alpha * 1.2, 0.3), o - 0.06), o]; };
  const TITLES = { bf: 'Virtual bright-field image', adf: 'Virtual annular dark-field image', vdf: 'Virtual dark-field image', ori: 'Crystal orientation map' };

  const scanBtn = btn('Scan', () => { if (!eng) return; if (scanning) scanning = false; else { if (acquired >= SG * SG) acquired = 0; scanning = true; } updStatus(); });
  const sA = slider({ label: 'Convergence semi-angle α', min: 0.04, max: 0.9, step: 0.01, value: 0.12, fmt: v => v.toFixed(2) + ' Å⁻¹ · ' + (v * ST.LAMBDA * 1000).toFixed(1) + ' mrad', on: v => { if (!eng) return; eng.alpha = v; eng.makeProbe(); restart(); } });
  const sD = slider({ label: 'Defocus Δf', min: -400, max: 400, step: 5, value: 0, fmt: v => (v > 0 ? '+' : '') + v + ' Å', on: v => { if (!eng) return; eng.df = v; eng.makeProbe(); restart(); } });
  const sG = slider({ label: 'Display gamma', min: 0.1, max: 1, step: 0.05, value: gamma, fmt: v => v.toFixed(2), on: v => { gamma = v; live(false); } });
  const sDet = seg({ label: 'Virtual detector', options: [['bf', 'Bright field'], ['adf', 'Annular DF'], ['vdf', 'Dark field (spot)'], ['ori', 'Orientation']], value: mode, on: v => { mode = v; if (v === 'vdf' && !vdf) vdf = brightestSpot(); remask(); } });
  const sView = seg({ label: 'Left panel shows', options: [['vimg', 'Virtual image'], ['sample', 'Atoms']], value: view, on: v => { view = v; titles(); dirty = true; } });

  root.append(
    h('div', { class: 'sim-scene' }, sceneCv),
    h('div', { class: 'sim-row' }, [
      h('div', { class: 'panel' }, [h('div', { class: 'panel-h' }, [aTitle, h('span', { text: 'drag to move probe' })]), aCv]),
      h('div', { class: 'panel' }, [h('div', { class: 'panel-h' }, [h('span', { text: 'Diffraction pattern at the probe' }), h('span', { text: 'click a spot → dark field' })]), bCv])]),
    h('div', { class: 'sim-bar' }, [scanBtn, status]),
    h('div', { class: 'ctrls' }, [sA, sD, sDet, sView, sG]));

  function titles() { aTitle.textContent = view === 'sample' ? 'Atoms, colored by grain' : TITLES[mode]; }
  function updStatus() {
    const n = SG * SG, pct = Math.round(acquired / n * 100);
    status.textContent = scanning ? 'Acquiring 4D dataset… ' + pct + '%' : acquired >= n ? '4D dataset: ' + SG + '×' + SG + ' probe positions × ' + ST.N + '×' + ST.N + '-pixel patterns' : acquired ? 'Paused at ' + pct + '%' : 'Press Scan to acquire';
    scanBtn.textContent = scanning ? 'Pause scan' : acquired >= n ? 'Rescan' : acquired ? 'Resume scan' : 'Scan';
    titles();
  }
  function computeMask() {
    const idx = [], dk = ST.dk, ar = adfRange();
    for (let by = 0; by < NB; by++) for (let bx = 0; bx < NB; bx++) {
      const kx = (2 * bx - 39.5) * dk, ky = (2 * by - 39.5) * dk, k = Math.hypot(kx, ky);
      const m = mode === 'bf' ? k <= eng.alpha : mode === 'adf' ? k >= ar[0] && k <= ar[1] : mode === 'vdf' && vdf ? Math.hypot(kx - vdf[0], ky - vdf[1]) <= vdfR : false;
      if (m) idx.push(by * NB + bx);
    }
    maskIdx = Int32Array.from(idx);
  }
  function updateVals(a, b) { for (let i = a; i < b; i++) { let s = 0; const o = i * NN; for (let j = 0; j < maskIdx.length; j++) s += cube[o + maskIdx[j]]; vals[i] = s; } }
  function paintV() {
    const d = vImg.data, n = SG * SG;
    if (mode === 'ori') { for (let i = 0; i < n; i++) { if (i < acquired) { const c = oriRGB(oriH[i], oriC[i]); d[i * 4] = c[0] * 255; d[i * 4 + 1] = c[1] * 255; d[i * 4 + 2] = c[2] * 255; d[i * 4 + 3] = 245; } else d[i * 4 + 3] = 0; } }
    else {
      let mn = Infinity, mx = -Infinity; for (let i = 0; i < acquired; i++) { if (vals[i] < mn) mn = vals[i]; if (vals[i] > mx) mx = vals[i]; }
      const sc = mx > mn ? 255 / (mx - mn) : 0;
      for (let i = 0; i < n; i++) { if (i < acquired) { const g = (vals[i] - mn) * sc; d[i * 4] = d[i * 4 + 1] = d[i * 4 + 2] = g; d[i * 4 + 3] = 245; } else d[i * 4 + 3] = 0; }
    }
    vctx.putImageData(vImg, 0, 0); composeStem(cctx, eng.S.tex, vCv, SG); dirty = true;
  }
  function live(recompute) { if (recompute !== false) eng.diffract(probe[0], probe[1]); paintDP(eng.I, dImg, gamma); dctx.putImageData(dImg, 0, 0); dirty = true; }
  function remask() { computeMask(); updateVals(0, acquired); paintV(); titles(); }
  function restart() { acquired = 0; scanning = true; computeMask(); live(); paintV(); updStatus(); }
  function brightestSpot() {
    const N = ST.N, c = N / 2, dk = ST.dk, r0 = (eng.alpha + 0.1) / dk, r1 = 1.1 / dk; let best = -1, bx = 0, by = 0;
    for (let y = 0; y < N; y++) for (let x = 0; x < N; x++) { const r = Math.hypot(x - c, y - c); if (r < r0 || r > r1) continue; const v = eng.I[y * N + x]; if (v > best) { best = v; bx = x; by = y; } }
    return best > 0 ? [(bx - c) * dk, (by - c) * dk] : [0.5, 0];
  }
  function setProbe(fx, fy) { probe = home = [Math.round(clamp(fx, ST.S0, ST.S1)), Math.round(clamp(fy, ST.S0, ST.S1))]; if (scanning) { scanning = false; updStatus(); } live(); }
  function stepScan() {
    const t0 = performance.now(), from = acquired, n = SG * SG;
    while (acquired < n && performance.now() - t0 < 10) {
      const p = scanPos(acquired, SG), I = eng.diffract(p[0], p[1]), o = acquired * NN, c0 = ST.N / 2 - NB, N = ST.N;
      for (let by = 0; by < NB; by++) { const y = c0 + 2 * by; for (let bx = 0; bx < NB; bx++) { const i = y * N + c0 + 2 * bx; cube[o + by * NB + bx] = I[i] + I[i + 1] + I[i + N] + I[i + N + 1]; } }
      const or = eng.orient(I); oriH[acquired] = or[0]; oriC[acquired] = or[1]; probe = p; acquired++;
    }
    updateVals(from, acquired);
    if (acquired >= n) { scanning = false; probe = home; live(); } else live(false); paintV(); updStatus();
  }
  function drawA() {
    const { ctx, w } = PA; PA.reset(); ctx.fillStyle = COL.panel; ctx.fillRect(0, 0, w, w);
    ctx.imageSmoothingEnabled = true; ctx.drawImage(view === 'vimg' ? comp : eng.S.tex, R0, R0, RS, RS, 0, 0, w, w);
    const sc = w / RS, X = (probe[0] - R0) * sc, Y = (probe[1] - R0) * sc;
    const dA = Math.max(0.51 / eng.alpha, 2 * Math.abs(eng.df) * eng.alpha * ST.LAMBDA), rp = Math.max(3, dA / ST.px * sc / 2);
    ctx.strokeStyle = COL.acc; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.arc(X, Y, rp, 0, TAU); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(X - rp - 6, Y); ctx.lineTo(X - rp - 1, Y); ctx.moveTo(X + rp + 1, Y); ctx.lineTo(X + rp + 6, Y); ctx.moveTo(X, Y - rp - 6); ctx.lineTo(X, Y - rp - 1); ctx.moveTo(X, Y + rp + 1); ctx.lineTo(X, Y + rp + 6); ctx.stroke();
    const nm = 10 / ST.px * sc; scaleBar(ctx, w - 14 - nm, w - 14, w - 14, '1 nm');
  }
  function drawB() {
    const { ctx, w } = PB, c = w / 2, s = w / (ST.N * ST.dk); PB.reset(); ctx.imageSmoothingEnabled = true; ctx.drawImage(dpCv, 0, 0, w, w);
    ctx.strokeStyle = 'rgba(120,235,100,.95)'; ctx.lineWidth = 1.5;
    const ring = (x, y, r) => { ctx.beginPath(); ctx.arc(x, y, r, 0, TAU); ctx.stroke(); };
    if (mode === 'bf') ring(c, c, eng.alpha * s);
    else if (mode === 'adf') { const ar = adfRange(); ctx.fillStyle = 'rgba(88,185,71,.12)'; ctx.beginPath(); ctx.arc(c, c, ar[1] * s, 0, TAU); ctx.arc(c, c, ar[0] * s, 0, TAU, true); ctx.fill(); ring(c, c, ar[0] * s); ring(c, c, ar[1] * s); }
    else if (mode === 'vdf' && vdf) ring(c + vdf[0] * s, c + vdf[1] * s, vdfR * s);
    else if (mode === 'ori') { ctx.setLineDash([4, 4]); ring(c, c, eng.oR[0] * s); ring(c, c, eng.oR[1] * s); ctx.setLineDash([]); }
    scaleBar(ctx, w - 14 - s, w - 14, w - 14, '1 Å⁻¹');
  }
  drag(PA, p => { const sc = RS / PA.w; if (eng) setProbe(R0 + p[0] * sc, R0 + p[1] * sc); });
  drag(SC, (p, ph) => { if (!eng) return; const q = SC.cam.unproject(p[0], p[1], 0); if (ph === 'down' && (Math.abs(q[0]) > ST.H || Math.abs(q[1]) > ST.H)) return; setProbe((q[0] + ST.H) / ST.px, (ST.H - q[1]) / ST.px); });
  aCv.addEventListener('keydown', e => { if (!eng) return; const m = { ArrowLeft: [-4, 0], ArrowRight: [4, 0], ArrowUp: [0, -4], ArrowDown: [0, 4] }[e.key]; if (m) { e.preventDefault(); setProbe(probe[0] + m[0], probe[1] + m[1]); } });
  bCv.addEventListener('click', e => {
    if (!eng) return; const p = pos(PB, e), N = ST.N, dk = ST.dk, s = PB.w / (N * dk), c = N / 2;
    let px = Math.round(c + (p[0] - PB.w / 2) / s / dk), py = Math.round(c + (p[1] - PB.w / 2) / s / dk), best = -1, bx = px, by = py;
    for (let y = py - 4; y <= py + 4; y++) for (let x = px - 4; x <= px + 4; x++) { if (x < 0 || y < 0 || x >= N || y >= N) continue; const v = eng.I[y * N + x]; if (v > best) { best = v; bx = x; by = y; } }
    vdf = [(bx - c) * dk, (by - c) * dk]; mode = 'vdf'; sDet.set('vdf'); remask();
  });
  updStatus();
  whenVisible(root, () => {
    if (!eng) { eng = new StemEngine(0.12); cube = new Float32Array(SG * SG * NN); oriH = new Float32Array(SG * SG); oriC = new Float32Array(SG * SG); vals = new Float32Array(SG * SG); computeMask(); live(); paintV(); scanning = true; updStatus(); }
    if (scanning) stepScan();
    if (dirty) { drawStemScene(SC, eng, comp, dpCv, probe, true); drawA(); drawB(); dirty = false; }
  });
}

/* ========================================================== FIB-SEM tomography */
const FB = { X: 6, Y: 4, Z: 3, PT: 0.16, YMAX: 3.2 }; // µm; PT = protective Pt cap
const FIB_RGB = { 1: [88, 185, 71], 2: [74, 158, 255], 3: [245, 184, 61] };
let FIB_VOL = null;
function fibVolume() {
  if (FIB_VOL) return FIB_VOL;
  const r = rng32(23), segs = [], pores = [], ags = [];
  // lithium filaments grow up from the Li anode (z = 0), wandering and branching
  function grow(x, y, z, dx, dy, dz, rad, depth) {
    for (let n = 0; n < 70 && z < FB.Z - 0.22; n++) {
      dx += gauss(r) * 0.32; dy += gauss(r) * 0.32; dz += 0.42; const m = Math.hypot(dx, dy, dz); dx /= m; dy /= m; dz /= m;
      const L = 0.12, nx = clamp(x + dx * L, 0.15, FB.X - 0.15), ny = clamp(y + dy * L, 0.12, FB.YMAX - 0.1), nz = z + dz * L;
      segs.push([x, y, z, nx, ny, nz, rad]); x = nx; y = ny; z = nz; rad = Math.max(0.03, rad * 0.986);
      if (depth < 3 && r() < 0.075) { const a = r() * TAU; grow(x, y, z, Math.cos(a), Math.sin(a), 0.3, rad * 0.82, depth + 1); }
      if (depth > 0 && r() < 0.035) break;
    }
  }
  [[1.3, 1.0], [3.3, 2.1], [4.7, 0.9], [2.4, 0.45], [5.2, 2.4]].forEach(p => grow(p[0], p[1], 0, 0, 0, 1, 0.09, 0));
  for (let i = 0; i < 22; i++) pores.push([0.3 + r() * 5.4, 0.15 + r() * 3.6, 0.3 + r() * 2.4, 0.06 + r() * 0.12]);
  for (let i = 0; i < 70; i++) {
    let x, y, z;
    if (i < 30) { const s = segs[(r() * segs.length) | 0], a = r() * TAU, b = r() * Math.PI, d = s[6] + 0.04 + r() * 0.08; x = s[3] + d * Math.sin(b) * Math.cos(a); y = s[4] + d * Math.sin(b) * Math.sin(a); z = s[5] + d * Math.cos(b); }
    else { x = 0.2 + r() * 5.6; y = 0.1 + r() * 3.8; z = 0.2 + r() * 2.6; }
    ags.push([x, y, clamp(z, 0.05, FB.Z - 0.05), 0.025 + r() * 0.03]);
  }
  return (FIB_VOL = { segs, pores, ags });
}
/* SEM image of the cross-section face at depth y, then threshold segmentation of that (noisy) image */
function fibSlice(y, W, opt) {
  const V = fibVolume(), Zt = FB.Z + FB.PT, ps = FB.X / W, H = Math.round(Zt / ps), lab = new Uint8Array(W * H), cap = Math.round(FB.PT / ps);
  lab.fill(4, 0, cap * W);
  const col = x => Math.floor(x / ps), row = z => Math.floor((Zt - z) / ps);
  for (const s of V.segs) {
    const r = s[6]; if (Math.min(s[1], s[4]) - r > y || Math.max(s[1], s[4]) + r < y) continue;
    const i0 = Math.max(0, col(Math.min(s[0], s[3]) - r)), i1 = Math.min(W - 1, col(Math.max(s[0], s[3]) + r)), j0 = Math.max(cap, row(Math.max(s[2], s[5]) + r)), j1 = Math.min(H - 1, row(Math.min(s[2], s[5]) - r));
    const ax = s[0], ay = s[1], az = s[2], bx = s[3] - ax, by = s[4] - ay, bz = s[5] - az, bb = bx * bx + by * by + bz * bz || 1e-9, r2 = r * r;
    for (let j = j0; j <= j1; j++) { const z = Zt - (j + 0.5) * ps; for (let i = i0; i <= i1; i++) {
      const px = (i + 0.5) * ps - ax, py = y - ay, pz = z - az; let t = (px * bx + py * by + pz * bz) / bb; t = t < 0 ? 0 : t > 1 ? 1 : t;
      const dx = px - t * bx, dy = py - t * by, dz = pz - t * bz; if (dx * dx + dy * dy + dz * dz <= r2 && lab[j * W + i] === 0) lab[j * W + i] = 1;
    } }
  }
  const sph = (arr, L, over) => { for (const p of arr) {
    const dy = p[1] - y; if (Math.abs(dy) >= p[3]) continue; const rr = Math.sqrt(p[3] * p[3] - dy * dy), r2 = rr * rr;
    const i0 = Math.max(0, col(p[0] - rr)), i1 = Math.min(W - 1, col(p[0] + rr)), j0 = Math.max(cap, row(p[2] + rr)), j1 = Math.min(H - 1, row(p[2] - rr));
    for (let j = j0; j <= j1; j++) { const z = Zt - (j + 0.5) * ps; for (let i = i0; i <= i1; i++) { const x = (i + 0.5) * ps; if ((x - p[0]) ** 2 + (z - p[2]) ** 2 <= r2) { const k = j * W + i; if (over || lab[k] <= 1) lab[k] = L; } } }
  } };
  sph(V.pores, 2, false); sph(V.ags, 3, true);
  const se = opt.se !== false, img = new Float32Array(W * H), rnd = rng32(opt.seed || 1), sig = 0.025 + 0.06 / Math.sqrt(opt.dwell || 1);
  for (let i = 0; i < W; i++) {
    let streak = 0; const cs = hash3(i, 11, 5) * 2 - 1, cs2 = hash3(i >> 1, 3, 9) * 2 - 1, x = (i + 0.5) * ps;
    for (let j = 0; j < H; j++) {
      const k = j * W + i, l = lab[k], z = Zt - (j + 0.5) * ps;
      const e = (i > 0 && lab[k - 1] !== l) || (i < W - 1 && lab[k + 1] !== l) || (j > 0 && lab[k - W] !== l) || (j < H - 1 && lab[k + W] !== l);
      let v;
      if (l === 4) v = (se ? 0.7 : 0.92) + 0.05 * (vnoise3(x * 9, y * 3, z * 9) - 0.5);
      else if (l === 0) v = (se ? 0.47 : 0.42) + 0.12 * (vnoise3(x / 0.3, y / 0.3, z / 0.3) - 0.5) + 0.05 * (vnoise3(x / 0.08, y / 0.08, z / 0.08) - 0.5);
      else if (l === 1) v = se ? 0.2 : 0.12; else if (l === 2) v = se ? 0.06 : 0.03; else v = se ? 0.9 : 1.0;
      if (e && l >= 1 && l <= 3) v += se ? [0, 0.2, 0.4, 0.05][l] : 0.05; // SE edge effect
      if (l === 2 || l === 3) streak = Math.min(0.13, streak + 0.015); else streak *= 0.99; // curtaining below voids/particles
      img[k] = v + streak * cs + 0.012 * cs2 + sig * gauss(rnd);
    }
  }
  // segmentation: 3x3 mean filter + global thresholds (what you'd actually do to the image)
  const seg = new Uint8Array(W * H), th = se ? [0.13, 0.33, 0.75] : [0.08, 0.26, 0.75];
  seg.fill(4, 0, cap * W);
  for (let j = cap; j < H; j++) for (let i = 0; i < W; i++) {
    let s = 0, n = 0; for (let b = -1; b <= 1; b++) for (let a = -1; a <= 1; a++) { const jj = j + b, ii = i + a; if (jj < cap || jj >= H || ii < 0 || ii >= W) continue; s += img[jj * W + ii]; n++; }
    const m = s / n; seg[j * W + i] = m < th[0] ? 2 : m < th[1] ? 1 : m > th[2] ? 3 : 0;
  }
  return { W, H, lab, seg, img, y };
}
function PtBuf() { this.a = new Float32Array(4 * 16384); this.n = 0; }
PtBuf.prototype.push = function (x, y, z, l) { if (this.n * 4 + 4 > this.a.length) { const b = new Float32Array(this.a.length * 2); b.set(this.a); this.a = b; } const o = this.n * 4; this.a[o] = x; this.a[o + 1] = y; this.a[o + 2] = z; this.a[o + 3] = l; this.n++; };
function fibPoints(sl, pb, stride) {
  const { W, H, seg, y } = sl, ps = FB.X / W, Zt = FB.Z + FB.PT;
  for (let j = 1; j < H - 1; j += stride) for (let i = 1; i < W - 1; i += stride) {
    const k = j * W + i, l = seg[k]; if (l === 0 || l === 4) continue;
    if (seg[k - 1] !== l || seg[k + 1] !== l || seg[k - W] !== l || seg[k + W] !== l) pb.push((i + 0.5) * ps, y, Zt - (j + 0.5) * ps, l);
  }
}
function fibImageTo(cv, sl, overlay) {
  const { W, H, img, seg } = sl; if (cv.width !== W || cv.height !== H) { cv.width = W; cv.height = H; }
  const ctx = cv.getContext('2d'), id = ctx.createImageData(W, H), d = id.data;
  for (let i = 0; i < W * H; i++) {
    const g = clamp(img[i], 0, 1) * 255; let r = g, gg = g, b = g;
    if (overlay) { const c = FIB_RGB[seg[i]]; if (c) { r = g * 0.45 + c[0] * 0.55; gg = g * 0.45 + c[1] * 0.55; b = g * 0.45 + c[2] * 0.55; } }
    d[i * 4] = r; d[i * 4 + 1] = gg; d[i * 4 + 2] = b; d[i * 4 + 3] = 255;
  }
  ctx.putImageData(id, 0, 0);
}
const FIB_ION_H = 2.6;
const fibESrc = ys => [FB.X / 2, ys - 3.4, FB.Z + FB.PT + 2.6];
function fibCamFit(cam, w, h, pad) {
  const Zt = FB.Z + FB.PT, pts = [];
  for (const x of [0, FB.X]) for (const y of [0, FB.Y]) for (const z of [-0.3, Zt]) pts.push([x, y, z]);
  for (const ys of [0, FB.YMAX]) { pts.push([0, ys, Zt + FIB_ION_H + 0.3], [FB.X, ys, Zt + FIB_ION_H + 0.3], fibESrc(ys)); }
  cam.fit(pts, w, h, pad);
}
function drawFibScene(V, s) {
  const { ctx, w, h } = V, cam = V.cam, P = (x, y, z) => cam.p(x, y, z), X = FB.X, Y = FB.Y, Zt = FB.Z + FB.PT, ys = s.ys;
  V.reset(); ctx.fillStyle = COL.panel; ctx.fillRect(0, 0, w, h); ctx.lineWidth = 1;
  const face = (pts, fill, stroke) => { poly(ctx, pts.map(q => P(q[0], q[1], q[2]))); ctx.fillStyle = fill; ctx.fill(); if (stroke) { ctx.strokeStyle = stroke; ctx.stroke(); } };
  // Li metal anode under the electrolyte
  face([[0, 0, 0], [X, 0, 0], [X, 0, -0.3], [0, 0, -0.3]], '#2c3a30', 'rgba(255,255,255,.06)');
  face([[X, 0, 0], [X, Y, 0], [X, Y, -0.3], [X, 0, -0.3]], '#223027', 'rgba(255,255,255,.06)');
  if (ys > 0) face([[0, 0, 0], [X, 0, 0], [X, ys, 0], [0, ys, 0]], '#33463a', 'rgba(255,255,255,.06)');
  // milled-away volume (ghost)
  if (ys > 0.01) { ctx.setLineDash([3, 4]); ctx.strokeStyle = 'rgba(255,255,255,.18)'; for (const e of [[[0, 0, Zt], [X, 0, Zt]], [[X, 0, Zt], [X, ys, Zt]], [[0, 0, Zt], [0, ys, Zt]], [[0, 0, 0], [0, 0, Zt]], [[X, 0, 0], [X, 0, Zt]]]) { const a = P(...e[0]), b = P(...e[1]); ctx.beginPath(); ctx.moveTo(a[0], a[1]); ctx.lineTo(b[0], b[1]); ctx.stroke(); } ctx.setLineDash([]); }
  // remaining block
  face([[X, ys, 0], [X, Y, 0], [X, Y, FB.Z], [X, ys, FB.Z]], '#1b2029');
  face([[X, ys, FB.Z], [X, Y, FB.Z], [X, Y, Zt], [X, ys, Zt]], '#5c6069');
  face([[0, ys, Zt], [X, ys, Zt], [X, Y, Zt], [0, Y, Zt]], '#777b84');
  const f00 = P(0, ys, Zt), f10 = P(X, ys, Zt), f01 = P(0, ys, 0), f11 = P(X, ys, 0);
  quad(V, s.prevCv && s.reveal < 1 ? s.prevCv : s.faceCv, f00, f10, f01);
  if (s.prevCv && s.reveal < 1 && s.reveal > 0) { const zc = Zt * (1 - s.reveal); ctx.save(); poly(ctx, [f00, f10, P(X, ys, zc), P(0, ys, zc)]); ctx.clip(); quad(V, s.faceCv, f00, f10, f01); ctx.restore(); const a = P(0, ys, zc), b = P(X, ys, zc); ctx.strokeStyle = 'rgba(140,240,120,.8)'; ctx.beginPath(); ctx.moveTo(a[0], a[1]); ctx.lineTo(b[0], b[1]); ctx.stroke(); }
  ctx.strokeStyle = 'rgba(255,255,255,.22)'; poly(ctx, [f00, f10, f11, f01]); ctx.stroke();
  // reconstruction so far, drawn in the volume it came from
  if (s.pts && s.pts.n) { const a = s.pts.a, n = s.pts.n, ev = s.every || 1; for (let l = 1; l <= 3; l++) { const c = FIB_RGB[l]; ctx.fillStyle = 'rgba(' + c + ',.85)'; ctx.beginPath(); for (let i = 0; i < n; i += ev) { const o = i * 4; if (a[o + 3] !== l) continue; const q = P(a[o], a[o + 1], a[o + 2]); ctx.rect(q[0] - 0.7, q[1] - 0.7, 1.4, 1.4); } ctx.fill(); } }
  // electron beam (images the face)
  const es = fibESrc(ys), tg = [P(...es)]; for (let i = 0; i < 20; i++) { const a = i / 20 * TAU; tg.push(P(X / 2 + 0.5 * Math.cos(a), ys, 1.5 + 0.5 * Math.sin(a))); }
  ctx.fillStyle = 'rgba(88,185,71,' + (0.1 + 0.25 * (s.ebeam || 0)) + ')'; poly(ctx, hull(tg)); ctx.fill();
  // ion beam (mills the next slice)
  const xi = s.ionX == null ? X / 2 : 0.15 + s.ionX * (X - 0.3), top = P(xi, ys, Zt + FIB_ION_H), b1 = P(xi - 0.07, ys, Zt), b2 = P(xi + 0.07, ys, Zt);
  ctx.fillStyle = s.ionX == null ? 'rgba(74,158,255,.22)' : 'rgba(74,158,255,.6)'; poly(ctx, [[top[0] - 2.5, top[1]], [top[0] + 2.5, top[1]], b2, b1]); ctx.fill();
  if (s.ionX != null) { const a = P(0.1, ys, Zt), b = P(xi, ys, Zt); ctx.strokeStyle = 'rgba(120,190,255,.9)'; ctx.lineWidth = 2; ctx.beginPath(); ctx.moveTo(a[0], a[1]); ctx.lineTo(b[0], b[1]); ctx.stroke(); ctx.lineWidth = 1; ctx.save(); ctx.shadowColor = '#6fb6ff'; ctx.shadowBlur = 12; ctx.fillStyle = '#d8ecff'; ctx.beginPath(); ctx.arc(b[0], b[1], 2.4, 0, TAU); ctx.fill(); ctx.restore(); }
  if (s.labels) {
    txt(ctx, 'ion beam · mills a slice', top[0] + 8, top[1] + 4, { color: '#8fc3ff' });
    const e2 = P(...es); txt(ctx, 'electron beam · images the face', e2[0] - 6, e2[1] - 10, { align: 'center', color: '#9be08c' });
    const an = P(X, Y * 0.6, -0.15); txt(ctx, 'Li metal', an[0] + 8, an[1] + 6);
    const el = P(X, Y * 0.75, FB.Z * 0.5); txt(ctx, 'solid electrolyte', el[0] + 8, el[1]);
  }
}
function drawRecon(V, pb, ysMax, cam) {
  const { ctx, w, h } = V, X = FB.X, Z = FB.Z, cx = X / 2, cy = FB.YMAX / 2, cz = Z / 2, R = Math.hypot(X, FB.YMAX, Z) / 2;
  V.reset(); ctx.fillStyle = COL.panel; ctx.fillRect(0, 0, w, h);
  cam.s = Math.min(w, h) * 0.47 / R; cam.cx = w / 2; cam.cy = h / 2 + 4;
  const P = (x, y, z) => cam.p(x - cx, y - cy, z - cz);
  const cs = []; for (const x of [0, X]) for (const y of [0, ysMax]) for (const z of [0, Z]) cs.push([x, y, z]);
  ctx.strokeStyle = 'rgba(255,255,255,.14)'; ctx.lineWidth = 1;
  for (let i = 0; i < 8; i++) for (let j = i + 1; j < 8; j++) { let d = 0; for (let c = 0; c < 3; c++) if (cs[i][c] !== cs[j][c]) d++; if (d !== 1) continue; const a = P(...cs[i]), b = P(...cs[j]); ctx.beginPath(); ctx.moveTo(a[0], a[1]); ctx.lineTo(b[0], b[1]); ctx.stroke(); }
  const NBk = 8, lists = []; for (let i = 0; i < NBk * 4; i++) lists.push([]);
  const a = pb.a, n = pb.n;
  for (let i = 0; i < n; i++) { const o = i * 4, x = a[o] - cx, y = a[o + 1] - cy, z = a[o + 2] - cz, d = cam.depth(x, y, z), b = clamp(Math.floor((d + R) / (2 * R) * NBk), 0, NBk - 1), q = cam.p(x, y, z); lists[b * 4 + a[o + 3]].push(q[0], q[1]); }
  for (let b = NBk - 1; b >= 0; b--) { const f = 1 - 0.6 * b / (NBk - 1); for (let l = 1; l <= 3; l++) { const L = lists[b * 4 + l]; if (!L.length) continue; const c = FIB_RGB[l]; ctx.fillStyle = 'rgb(' + (c[0] * f | 0) + ',' + (c[1] * f | 0) + ',' + (c[2] * f | 0) + ')'; ctx.beginPath(); for (let i = 0; i < L.length; i += 2) ctx.rect(L[i] - 0.9, L[i + 1] - 0.9, 1.8, 1.8); ctx.fill(); } }
  const ax = P(0, 0, 0), ex = P(1, 0, 0), ey = P(0, 1, 0), ez = P(0, 0, 1);
  [[ex, 'x'], [ey, 'y (mill)'], [ez, 'z']].forEach(([e, t]) => { ctx.strokeStyle = 'rgba(255,255,255,.5)'; ctx.beginPath(); ctx.moveTo(ax[0], ax[1]); ctx.lineTo(e[0], e[1]); ctx.stroke(); txt(ctx, t, e[0] + (e[0] - ax[0]) * 0.25, e[1] + (e[1] - ax[1]) * 0.25, { size: 10, align: 'center' }); });
  txt(ctx, n ? sci(n) + ' surface points' : 'mill to start reconstructing', 10, h - 12, { size: 10 });
}

function fibPillar(cv) {
  const V = makeCanvas(cv, 360, 360); V.cam = new Cam(-32, 26); fibCamFit(V.cam, 360, 360, 12);
  const W = 150, t = 0.1, n = Math.floor(FB.YMAX / t + 1e-6); let s = null;
  function reset() { s.k = 0; s.pts = new PtBuf(); s.phase = 'mill'; s.tp = 0; s.reveal = 1; s.hold = 0; const sl = fibSlice(0, W, { seed: 3 }); fibPoints(sl, s.pts, 1); fibImageTo(s.face, sl, false); }
  function advance() { s.pctx.drawImage(s.face, 0, 0); s.k++; const sl = fibSlice(s.k * t, W, { seed: s.k * 31 + 3 }); fibPoints(sl, s.pts, 1); fibImageTo(s.face, sl, false); }
  whenVisible(cv, dt => {
    if (!s) { const face = offscreen(W, 80), prev = offscreen(W, 80); s = { face, prev, pctx: null }; reset(); prev.width = face.width; prev.height = face.height; s.pctx = prev.getContext('2d'); if (REDUCED) { for (let i = 0; i < 18; i++) advance(); } }
    if (REDUCED) { if (!s.done) { drawFibScene(V, { ys: s.k * t, faceCv: s.face, pts: s.pts, every: 2 }); s.done = true; } return; }
    let ionX = null, eb = 0;
    if (s.hold > 0) { s.hold -= dt; if (s.hold <= 0) reset(); }
    else { s.tp += dt; if (s.phase === 'mill') { ionX = Math.min(1, s.tp / 0.2); if (s.tp >= 0.2) { advance(); s.phase = 'img'; s.tp = 0; } } else { s.reveal = Math.min(1, s.tp / 0.14); eb = 1; if (s.tp >= 0.14) { s.reveal = 1; s.tp = 0; s.phase = 'mill'; if (s.k >= n) s.hold = 2.2; } } }
    drawFibScene(V, { ys: s.k * t, faceCv: s.face, prevCv: s.prev, reveal: s.phase === 'img' ? s.reveal : 1, ionX, ebeam: eb, pts: s.pts, every: 2 });
  });
}

function fibWidget(root) {
  const W = 240, Zt = FB.Z + FB.PT;
  const sceneCv = h('canvas', { class: 'scene', role: 'img', 'aria-label': '3D view of FIB-SEM slice-and-view: the ion beam mills slices from a solid-electrolyte block while the electron beam images each new face. Segmented features appear in the milled-away volume.' });
  const sCv = h('canvas', { 'aria-label': 'SEM image of the current cross-section face' }), rCv = h('canvas', { 'aria-label': '3D reconstruction. Drag to rotate.' });
  const SC = makeCanvas(sceneCv, 760, 420), PS = makeCanvas(sCv, 360, Math.round(360 * Zt / FB.X)), PR = makeCanvas(rCv, 360, 300);
  SC.cam = new Cam(-32, 26); fibCamFit(SC.cam, 760, 420, 18);
  const rcam = new Cam(-35, 24); let rDrag = null, idle = 0;
  let thick = 0.05, nSl, k, playing = false, phase = 'idle', tp = 0, reveal = 1, se = true, dwell = 1, overlay = false, review = -1, pts, cur, dirty = true, started = false;
  const face = offscreen(W, 10), prev = offscreen(W, 10), viewCv = offscreen(W, 10);
  const status = h('span', { class: 'hint', 'aria-live': 'polite' });
  const playBtn = btn('Mill ▶', () => { if (k >= nSl) reset(); playing = !playing; if (playing) { review = -1; sRev.set(0); if (phase === 'idle') { phase = 'mill'; tp = 0; } } upd(); });
  const finBtn = btn('Finish', () => { while (k < nSl) advance(); playing = false; phase = 'idle'; reveal = 1; upd(); });
  const resetBtn = btn('Reset', () => reset(), 'ghost');
  const sT = slider({ label: 'Slice thickness', min: 20, max: 150, step: 10, value: 50, fmt: v => v + ' nm', on: v => { thick = v / 1000; reset(); } });
  const sW = slider({ label: 'Dwell time per pixel', min: 0.1, max: 3, step: 0.1, value: 1, fmt: v => v.toFixed(1) + ' µs', on: v => { dwell = v; } });
  const sDet = seg({ label: 'SEM detector', options: [['se', 'Secondary (SE)'], ['bse', 'Backscattered (BSE)']], value: 'se', on: v => { se = v === 'se'; } });
  const sOv = seg({ label: 'Slice view', options: [['raw', 'Raw image'], ['seg', 'Segmented']], value: 'raw', on: v => { overlay = v === 'seg'; refreshFace(); } });
  const sRev = slider({ label: 'Review an earlier slice', min: 0, max: 0, step: 1, value: 0, fmt: v => review < 0 ? 'latest' : '#' + v, on: v => { review = v; playing = false; upd(); refreshView(); } });
  root.append(
    h('div', { class: 'sim-scene' }, sceneCv),
    h('div', { class: 'sim-row' }, [
      h('div', { class: 'panel' }, [h('div', { class: 'panel-h' }, [h('span', { text: 'SEM image of the face' }), h('span', { class: 'sl-lab' })]), sCv]),
      h('div', { class: 'panel' }, [h('div', { class: 'panel-h' }, [h('span', { text: '3D reconstruction' }), h('span', { text: 'drag to rotate' })]), rCv,
        h('div', { class: 'legend' }, [h('span', { class: 'lg', style: '--c:#58b947', text: 'Li filament' }), h('span', { class: 'lg', style: '--c:#4a9eff', text: 'void' }), h('span', { class: 'lg', style: '--c:#f5b83d', text: 'Ag-rich particle' })])])]),
    h('div', { class: 'sim-bar' }, [playBtn, finBtn, resetBtn, status]),
    h('div', { class: 'ctrls' }, [sT, sW, sDet, sOv, sRev]));
  const slLab = root.querySelector('.sl-lab');
  const opts = i => ({ se, dwell, seed: i * 7919 + 13 });
  function refreshFace() { if (cur) fibImageTo(face, cur, overlay); refreshView(); }
  function refreshView() {
    const i = review < 0 ? k : Math.min(review, k), sl = i === k && cur ? cur : fibSlice(i * thick, W, opts(i));
    fibImageTo(viewCv, sl, overlay); slLab.textContent = 'slice ' + i + ' of ' + nSl + ' · depth ' + (i * thick).toFixed(2) + ' µm'; dirty = true;
  }
  function reset() {
    nSl = Math.floor(FB.YMAX / thick + 1e-6); k = 0; pts = new PtBuf(); playing = false; phase = 'idle'; tp = 0; reveal = 1; review = -1;
    cur = fibSlice(0, W, opts(0)); fibPoints(cur, pts, 1); fibImageTo(face, cur, overlay); prev.width = face.width; prev.height = face.height;
    const r = sRev.querySelector('input'); r.max = 0; r.value = 0; sRev.set(0); refreshView(); upd();
  }
  function advance() {
    prev.getContext('2d').drawImage(face, 0, 0); k++; cur = fibSlice(k * thick, W, opts(k)); fibPoints(cur, pts, 1); fibImageTo(face, cur, overlay);
    sRev.querySelector('input').max = k; if (review < 0) refreshView();
  }
  function upd() {
    playBtn.textContent = playing ? 'Pause ❚❚' : k >= nSl ? 'Restart ▶' : k ? 'Resume ▶' : 'Mill ▶';
    status.textContent = k >= nSl ? 'Done: ' + (nSl + 1) + ' slices, ' + (FB.YMAX).toFixed(1) + ' µm milled at ' + Math.round(thick * 1000) + ' nm per slice' : 'Slice ' + k + ' of ' + nSl;
    dirty = true;
  }
  drag(PR, (p, ph) => { if (ph === 'down') rDrag = { p, yaw: rcam.yaw, el: rcam.elev }; else if (ph === 'move' && rDrag) { rcam.set(rDrag.yaw - (p[0] - rDrag.p[0]) * 0.5, clamp(rDrag.el + (p[1] - rDrag.p[1]) * 0.4, -10, 85)); dirty = true; } else rDrag = null; idle = 0; });
  whenVisible(root, dt => {
    if (!started) { started = true; reset(); }
    if (playing) {
      tp += dt;
      if (phase === 'mill' && tp >= 0.22) { advance(); phase = 'img'; tp = 0; upd(); }
      else if (phase === 'img') { reveal = Math.min(1, tp / 0.14); if (tp >= 0.14) { reveal = 1; tp = 0; phase = 'mill'; if (k >= nSl) { playing = false; phase = 'idle'; upd(); } } }
      dirty = true;
    }
    if (!rDrag && !REDUCED) { idle += dt; if (idle > 1.5) { rcam.set(rcam.yaw + dt * 10, rcam.elev); dirty = true; } }
    if (!dirty) return;
    drawFibScene(SC, { ys: k * thick, faceCv: face, prevCv: prev, reveal: phase === 'img' ? reveal : 1, ionX: playing && phase === 'mill' ? Math.min(1, tp / 0.22) : null, ebeam: phase === 'img' ? 1 : 0, pts, every: 3, labels: true });
    const { ctx, w, h: hh } = PS; PS.reset(); ctx.imageSmoothingEnabled = true; ctx.drawImage(viewCv, 0, 0, w, hh);
    scaleBar(ctx, w - 14 - w / FB.X, w - 14, hh - 12, '1 µm');
    drawRecon(PR, pts, k * thick, rcam); dirty = false;
  });
}

/* ================================================================== EELS */
/* sample: Li metal | interphase (Li2O / Li2CO3) | garnet LLZO (Li7La3Zr2O12), 120 x 60 nm */
const EX = 120, EY = 60, EZ = { Li: 3, C: 6, O: 8, Zr: 40, La: 57 };
const CARB = [[44, 12, 3.5], [47, 38, 4], [42, 52, 3]];
const MAT = [{ Ep: 7.1, G: 2.2, b: 0.55 }, { Ep: 21.5, G: 7, b: 0.9 }, { Ep: 23.5, G: 8.5, b: 1.2 }];
function drudeRaw(E, Ep, G) { if (E <= 0) return 0; const d = E * E - Ep * Ep; return E * G * Ep * Ep / (d * d + E * E * G * G) / (1 + Math.pow(E / (3 * Ep), 4)); }
MAT.forEach(m => { let s = 0; const d = m.Ep / 300; for (let E = d / 2; E < m.Ep * 25; E += d) s += drudeRaw(E, m.Ep, m.G) * d; m.norm = s; });
const Bstd = E => 0.02 * Math.pow(E / 20, -3);
const EDGE = { Li: { E0: 55, R: 3, lab: 'Li K' }, Zr: { E0: 181, R: 2, lab: 'Zr M₄,₅' }, C: { E0: 284, R: 9, lab: 'C K' }, O: { E0: 532, R: 6, lab: 'O K' }, La: { E0: 832, R: 16, lab: 'La M₄,₅' } };
const ELS = ['Li', 'Zr', 'C', 'O', 'La'];
ELS.forEach(k => { EDGE[k].amp = EDGE[k].R * Bstd(EDGE[k].E0); });
const FITW = { Li: [44, 53, 55, 70], Zr: [150, 176, 181, 230], C: [258, 281, 284, 314], O: [500, 528, 532, 562], La: [795, 826, 828, 858] };
const MAP_EL = [['Li', [88, 185, 71], 0.25], ['C', [245, 184, 61], 0.5], ['O', [255, 107, 90], 0.5], ['La', [74, 158, 255], 0.25]];
const WIN = { ll: { name: 'Low loss', e0: -3, e1: 80, d: 0.1, log: true, edges: ['Li'] }, li: { name: 'Li K', e0: 38, e1: 110, d: 0.2, edges: ['Li'] }, c: { name: 'C K', e0: 245, e1: 345, d: 0.25, edges: ['C'] }, o: { name: 'O K', e0: 490, e1: 600, d: 0.25, edges: ['O'] }, la: { name: 'La M₄,₅', e0: 785, e1: 885, d: 0.2, edges: ['La'] } };
function eelsComp(x, y) {
  const b1 = 36 + 3.5 * Math.sin(y / 7.5) + 1.8 * Math.sin(y / 3.1 + 1), b2 = 58 + 4 * Math.sin(y / 10 + 2) + 2 * Math.sin(y / 4.3);
  const s1 = 1 / (1 + Math.exp(-(x - b1) / 1.2)), s2 = 1 / (1 + Math.exp(-(x - b2) / 1.5));
  const fLi = 1 - s1, fLL = s2, fIn = s1 * (1 - s2);
  let carb = 0; for (const c of CARB) carb = Math.max(carb, 1 / (1 + Math.exp((Math.hypot(x - c[0], y - c[1]) - c[2]) / 0.8))); carb *= fIn;
  const el = { Li: fLi + fIn * 0.5 + fLL * 0.29, O: fIn * 0.38 + fLL * 0.5, C: fIn * (0.04 + 0.22 * carb), La: fLL * 0.125, Zr: fLL * 0.083 };
  const f = [fLi, fIn, fLL];
  return { f, el, b: f[0] * MAT[0].b + f[1] * MAT[1].b + f[2] * MAT[2].b, Ep: f[0] * MAT[0].Ep + f[1] * MAT[1].Ep + f[2] * MAT[2].Ep };
}
const smoothStep = (E, E0, s) => 0.5 * (1 + erf((E - E0) / (s * 1.41421)));
const gpk = (E, mu, w, s) => { const sg = Math.sqrt(w * w + s * s); return Math.exp(-((E - mu) * (E - mu)) / (2 * sg * sg)); };
function edgeShape(el, E, s) {
  switch (el) {
    case 'Li': return smoothStep(E, 55, s + 0.4) * Math.pow(Math.max(E, 55) / 55, -2.6) + 0.5 * gpk(E, 62, 3, s);
    case 'Zr': return smoothStep(E, 181, s + 1.5) * 0.35 * Math.pow(Math.max(E, 181) / 181, -2) + 0.9 * gpk(E, 215, 18, 0) / (1 + Math.exp(-(E - 185) / 4));
    case 'C': return smoothStep(E, 284, s + 0.3) * 0.55 * Math.pow(Math.max(E, 284) / 284, -2) + 1.2 * gpk(E, 290.3, 0.7, s) + 0.6 * gpk(E, 301, 4, s);
    case 'O': return smoothStep(E, 532, s + 0.4) * 0.6 * Math.pow(Math.max(E, 532) / 532, -1.8) + 0.7 * gpk(E, 534.5, 1.1, s) + 0.9 * gpk(E, 541, 4.5, s);
    default: return smoothStep(E, 830, s + 0.5) * 0.1 * Math.pow(Math.max(E, 830) / 830, -1.5) + gpk(E, 832, 1.0, s) + 1.3 * gpk(E, 849, 1.1, s);
  }
}
/* single-scattering + Poisson multiple scattering (plasmons), power-law background, core-loss edges */
function eelsModel(E, c, P) {
  const t = P.t, s = Math.max(0.05, P.res) / 2.355; let v = 0;
  if (E < 200) {
    v += Math.exp(-t) * Math.exp(-E * E / (2 * s * s)) / (s * 2.5066);
    let tn = 1, fa = 1;
    for (let n = 1; n <= 4; n++) { tn *= t; fa *= n; const wn = Math.exp(-t) * tn / fa; for (let m = 0; m < 3; m++) { const f = c.f[m]; if (f < 1e-3) continue; const M = MAT[m]; let sh; if (n === 1) sh = drudeRaw(E, M.Ep, M.G) / M.norm; else { const sg = Math.sqrt(n) * M.G / 2.355 * 1.15; sh = Math.exp(-(E - n * M.Ep) * (E - n * M.Ep) / (2 * sg * sg)) / (sg * 2.5066); } v += wn * f * sh; } }
  }
  if (E > 2) {
    let core = Bstd(E) * c.b * (1 + 0.6 * t) * smoothStep(E, 28, 6);
    for (const el of ELS) { const cc = c.el[el], E0 = EDGE[el].E0; if (cc < 1e-4 || E < E0 - 25) continue; core += EDGE[el].amp * cc * (edgeShape(el, E, s) + (E > E0 + c.Ep - 25 ? 0.8 * t * edgeShape(el, E - c.Ep, s + c.Ep * 0.12) : 0)); }
    v += t * Math.exp(-t) * core;
  }
  return v;
}
function noisy(l) { if (l <= 0) return 0; if (l < 30) { const L = Math.exp(-l); let k = 0, p = 1; do { k++; p *= Math.random(); } while (p > L); return k - 1; } return Math.max(0, l + Math.sqrt(l) * gauss(Math.random)); }
function plfit(E, Y, a, b) {
  let n = 0, sx = 0, sy = 0, sxx = 0, sxy = 0;
  for (let k = 0; k < E.length; k++) if (E[k] >= a && E[k] <= b && Y[k] > 0) { const x = Math.log(E[k]), y = Math.log(Y[k]); n++; sx += x; sy += y; sxx += x * x; sxy += x * y; }
  if (n < 3) return null; const m = (n * sxy - sx * sy) / (n * sxx - sx * sx); return { A: Math.exp((sy - m * sx) / n), r: -m };
}
function edgeSignal(el, E, Y) { const w = FITW[el], f = plfit(E, Y, w[0], w[1]); if (!f) return { sig: 0 }; let s = 0; for (let k = 0; k < E.length; k++) if (E[k] >= w[2] && E[k] <= w[3]) s += Y[k] - f.A * Math.pow(E[k], -f.r); return { sig: s, A: f.A, r: f.r }; }
function eelsAcquire(c, P, dose, e0, e1, d) { const n = Math.round((e1 - e0) / d) + 1, E = new Float64Array(n), Y = new Float64Array(n); for (let k = 0; k < n; k++) { E[k] = e0 + k * d; Y[k] = noisy(dose * d * eelsModel(E[k], c, P)); } return { E, Y }; }
function eelsTexture() {
  const W = 240, H = 120, cv = offscreen(W, H), ctx = cv.getContext('2d'), id = ctx.createImageData(W, H), d = id.data, r = rng32(5);
  for (let j = 0; j < H; j++) for (let i = 0; i < W; i++) {
    const x = (i + 0.5) * EX / W, y = (j + 0.5) * EY / H, c = eelsComp(x, y); let s = 0; for (const el in c.el) s += c.el[el] * Math.pow(EZ[el], 1.7);
    let v = 0.8 * Math.pow(s / 185, 0.42) * (1 + 0.22 * c.f[2] * (vnoise3(x / 6, y / 6, 1.5) - 0.5)) + 0.035 * gauss(r); v = clamp(v, 0, 1) * 255;
    const k = (j * W + i) * 4; d[k] = v * 0.96; d[k + 1] = v; d[k + 2] = v * 1.04; d[k + 3] = 255;
  }
  ctx.putImageData(id, 0, 0); return cv;
}
function paintOverview(cv, c, P) {
  const ctx = cv.getContext('2d'), w = cv.width, h = cv.height, e1 = 900, lmin = -7.4, lmax = 0.6, pad = 14;
  ctx.fillStyle = '#0b0e13'; ctx.fillRect(0, 0, w, h);
  ctx.strokeStyle = 'rgba(255,255,255,.08)'; ctx.lineWidth = 1; for (let i = 1; i < 4; i++) { const y = pad + (h - 2 * pad) * i / 4; ctx.beginPath(); ctx.moveTo(pad, y); ctx.lineTo(w - pad, y); ctx.stroke(); }
  const X = E => pad + (w - 2 * pad) * E / e1, Yv = v => h - pad - (h - 2 * pad) * (clamp(Math.log10(Math.max(v, 1e-9)), lmin, lmax) - lmin) / (lmax - lmin);
  const pts = []; for (let i = 0; i <= w - 2 * pad; i++) { const Ea = i / (w - 2 * pad) * e1, Eb = (i + 1) / (w - 2 * pad) * e1; let m = 0; const sub = Ea < 120 ? 24 : 2; for (let q = 0; q < sub; q++) m = Math.max(m, eelsModel(Ea + (Eb - Ea) * q / sub, c, P)); pts.push([pad + i, Yv(m)]); }
  ctx.beginPath(); ctx.moveTo(pad, h - pad); pts.forEach(p => ctx.lineTo(p[0], p[1])); ctx.lineTo(w - pad, h - pad); ctx.closePath(); ctx.fillStyle = 'rgba(88,185,71,.16)'; ctx.fill();
  ctx.beginPath(); pts.forEach((p, i) => (i ? ctx.lineTo(p[0], p[1]) : ctx.moveTo(p[0], p[1]))); ctx.strokeStyle = '#7fdc6c'; ctx.lineWidth = 1.6; ctx.stroke();
  ctx.font = '600 11px ' + FONT; ctx.textAlign = 'center';
  for (const el of ['Li', 'C', 'O', 'La']) { const E0 = EDGE[el].E0, on = c.el[el] > 0.02; ctx.fillStyle = on ? '#f0ede8' : 'rgba(240,237,232,.3)'; ctx.fillText(el, X(E0) + 6, Yv(eelsModel(E0 + 6, c, P)) - 9); }
  ctx.fillStyle = 'rgba(240,237,232,.55)'; ctx.font = '500 10px ' + FONT; ctx.textAlign = 'right'; ctx.fillText('energy loss →', w - pad, h - 4);
}
function eelsCamFit(cam, w, h, pad) { const pts = []; for (const x of [-EX / 2, EX / 2]) for (const y of [-EY / 2, EY / 2]) pts.push([x, y, 0], [x, y, -3]); pts.push([0, 0, 58], [-12, -11, -42], [78, -42, -68], [78, 42, -10], [78, -42, -10], [78, 42, -68]); cam.fit(pts, w, h, pad); }
function drawEelsScene(V, s) {
  const { ctx, w, h } = V, cam = V.cam, P = (x, y, z) => cam.p(x, y, z), wx = s.probe[0] - EX / 2, wy = EY / 2 - s.probe[1];
  V.reset(); ctx.fillStyle = COL.panel; ctx.fillRect(0, 0, w, h); ctx.lineWidth = 1;
  const face = (pts, fill, stroke) => { poly(ctx, pts.map(q => P(q[0], q[1], q[2]))); ctx.fillStyle = fill; ctx.fill(); if (stroke) { ctx.strokeStyle = stroke; ctx.stroke(); } };
  // spectrometer: prism, dispersed fan, detector with the spectrum
  const D = 78, ex = P(20, 0, -38), lo = P(D, -40, -38), hi = P(D, 40, -38);
  ctx.fillStyle = 'rgba(88,185,71,.10)'; cone(ctx, cam, [wx, wy, -3, 0.8], [0, 0, -34, 1.5]);
  face([[-9, -11, -34], [20, -11, -34], [20, 11, -34], [-9, 11, -34]], '#2a303c', 'rgba(255,255,255,.12)');
  face([[20, -11, -34], [20, 11, -34], [20, 11, -42], [20, -11, -42]], '#1f2430', 'rgba(255,255,255,.1)');
  face([[-9, -11, -34], [20, -11, -34], [20, -11, -42], [-9, -11, -42]], '#262b37', 'rgba(255,255,255,.1)');
  const g = ctx.createLinearGradient(lo[0], lo[1], hi[0], hi[1]); g.addColorStop(0, 'rgba(120,235,100,.38)'); g.addColorStop(0.5, 'rgba(90,200,200,.22)'); g.addColorStop(1, 'rgba(74,158,255,.18)');
  ctx.fillStyle = g; poly(ctx, [P(20, -1, -38), P(20, 1, -38), hi, lo]); ctx.fill();
  for (let i = 0; i <= 6; i++) { const q = P(D, -40 + i * 80 / 6, -38); ctx.strokeStyle = 'rgba(160,230,200,' + (0.35 - i * 0.03) + ')'; ctx.beginPath(); ctx.moveTo(ex[0], ex[1]); ctx.lineTo(q[0], q[1]); ctx.stroke(); }
  const s00 = P(D, -42, -10), s10 = P(D, 42, -10), s01 = P(D, -42, -68), s11 = P(D, 42, -68);
  quad(V, s.specCv, s00, s10, s01); ctx.strokeStyle = 'rgba(255,255,255,.2)'; poly(ctx, [s00, s10, s11, s01]); ctx.stroke();
  // lamella
  const H2 = EY / 2, W2 = EX / 2, T = [[-W2, H2], [W2, H2], [W2, -H2], [-W2, -H2]], sides = [];
  for (let i = 0; i < 4; i++) { const a = T[i], b = T[(i + 1) % 4]; sides.push({ q: [P(a[0], a[1], 0), P(b[0], b[1], 0), P(b[0], b[1], -3), P(a[0], a[1], -3)], d: cam.depth((a[0] + b[0]) / 2, (a[1] + b[1]) / 2, -1.5), i }); }
  sides.sort((u, v) => v.d - u.d); for (const sd of sides) { ctx.fillStyle = sd.i % 2 ? '#141923' : '#1b2230'; poly(ctx, sd.q); ctx.fill(); }
  const t00 = P(-W2, H2, 0), t10 = P(W2, H2, 0), t01 = P(-W2, -H2, 0), t11 = P(W2, -H2, 0);
  quad(V, s.tex, t00, t10, t01); ctx.strokeStyle = 'rgba(255,255,255,.2)'; poly(ctx, [t00, t10, t11, t01]); ctx.stroke();
  // probe
  ctx.fillStyle = 'rgba(88,185,71,.2)'; ctx.strokeStyle = 'rgba(130,225,110,.5)'; cone(ctx, cam, [wx, wy, 52, 9], [wx, wy, 0, 0.6], true);
  ctx.strokeStyle = 'rgba(170,245,150,.75)'; poly(ctx, circ3(wx, wy, 52, 9, 40).map(q => P(q[0], q[1], q[2]))); ctx.stroke();
  const sp = P(wx, wy, 0); ctx.save(); ctx.shadowColor = '#7dff6a'; ctx.shadowBlur = 14; ctx.fillStyle = '#eaffe4'; ctx.beginPath(); ctx.arc(sp[0], sp[1], 2.6, 0, TAU); ctx.fill(); ctx.restore();
  if (s.labels) {
    const pt = P(wx, wy, 52); txt(ctx, 'electron probe', pt[0], pt[1] - 9 * cam.s * cam.sE - 10, { align: 'center', color: COL.t1 });
    const pr = P(-9, -11, -42); txt(ctx, 'magnetic prism', pr[0] - 6, pr[1] + 10, { align: 'right' });
    txt(ctx, 'energy-loss spectrum', s01[0], s01[1] + 14, { align: 'left' });
    const lm = P(-W2, -H2, -3); txt(ctx, 'sample', lm[0] - 6, lm[1] + 4, { align: 'right' });
  }
}
function eelsPillar(cv) {
  const V = makeCanvas(cv, 360, 360); V.cam = new Cam(-28, 26); eelsCamFit(V.cam, 360, 360, 10);
  let s = null; const P = { t: 0.45, res: 0.5 };
  whenVisible(cv, dt => {
    if (!s) s = { tex: eelsTexture(), spec: offscreen(320, 170), t: 1.2 };
    if (REDUCED && s.done) return;
    s.t += dt; const px = 60 + 46 * Math.sin(s.t * 0.4), probe = [px, 32];
    paintOverview(s.spec, eelsComp(px, 32), P); drawEelsScene(V, { probe, tex: s.tex, specCv: s.spec }); s.done = true;
  });
}
function eelsWidget(root) {
  const MX = 48, MY = 24;
  const sceneCv = h('canvas', { class: 'scene', role: 'img', 'aria-label': '3D view of an EELS experiment: the probe passes through a lamella, then a magnetic prism disperses the electrons by energy onto a detector showing the spectrum.' });
  const aCv = h('canvas', { tabindex: '0', 'aria-label': 'Sample image. Drag, or use the arrow keys, to move the probe.' }), bCv = h('canvas', { 'aria-label': 'Energy-loss spectrum' });
  const SC = makeCanvas(sceneCv, 760, 380), PA = makeCanvas(aCv, 360, 180), PB = makeCanvas(bCv, 420, 270);
  SC.cam = new Cam(-28, 24); eelsCamFit(SC.cam, 760, 380, 18);
  let probe = [47, 30], P = { t: 0.4, res: 0.3 }, dose = 1e10, win = 'o', logY = false, view = 'haadf', spec = null, dirty = true, started = false, acq = 0, mapping = false;
  const show = { Li: true, C: true, O: true, La: true }, maps = {}; MAP_EL.forEach(m => { maps[m[0]] = new Float32Array(MX * MY); });
  let haadf = null; const overview = offscreen(320, 170), comp = offscreen(240, 120), cctx = comp.getContext('2d'), mapCv = offscreen(MX, MY), mctx = mapCv.getContext('2d'), mImg = mctx.createImageData(MX, MY);
  const readout = h('div', { class: 'readout', 'aria-live': 'polite' }), status = h('span', { class: 'hint', 'aria-live': 'polite' });
  const sWin = seg({ label: 'Energy window', options: Object.keys(WIN).map(k => [k, WIN[k].name]), value: win, on: v => { win = v; logY = !!WIN[v].log; sLog.set(logY ? 'log' : 'lin'); acquire(); } });
  const sLog = seg({ label: 'Intensity scale', options: [['lin', 'Linear'], ['log', 'Log']], value: 'lin', on: v => { logY = v === 'log'; dirty = true; } });
  const sT = slider({ label: 'Thickness t/λ', min: 0.1, max: 2, step: 0.05, value: P.t, fmt: v => v.toFixed(2), on: v => { P.t = v; changed(); } });
  const sR = slider({ label: 'Energy resolution', min: 0.1, max: 2, step: 0.05, value: P.res, fmt: v => v.toFixed(2) + ' eV', on: v => { P.res = v; changed(); } });
  const sE = slider({ label: 'Exposure (electrons)', min: 0, max: 2, step: 0.1, value: 1, fmt: v => sci(1e9 * Math.pow(10, v)), on: v => { dose = 1e9 * Math.pow(10, v); changed(); } });
  const mapBtn = btn('Acquire spectrum image', () => { if (mapping) { mapping = false; } else { if (acq >= MX * MY) acq = 0; mapping = true; view = 'map'; sView.set('map'); } upd(); });
  const sView = seg({ label: 'Sample view', options: [['haadf', 'HAADF image'], ['map', 'EELS map']], value: view, on: v => { view = v; compose(); } });
  const chips = h('div', { class: 'chips' }, MAP_EL.map(([el, c]) => h('button', { type: 'button', class: 'chip', 'aria-pressed': 'true', style: '--c:rgb(' + c + ')', text: EDGE[el].lab, onclick: e => { show[el] = !show[el]; e.currentTarget.setAttribute('aria-pressed', String(show[el])); paintMap(); } })));
  root.append(
    h('div', { class: 'sim-scene' }, sceneCv),
    h('div', { class: 'sim-row eels' }, [
      h('div', {}, [h('div', { class: 'panel' }, [h('div', { class: 'panel-h' }, [h('span', { text: 'Li metal · interphase · LLZO' }), h('span', { text: 'drag to move probe' })]), aCv]), chips, h('div', { class: 'sim-bar' }, [mapBtn, status])]),
      h('div', {}, [h('div', { class: 'panel' }, [h('div', { class: 'panel-h' }, [h('span', { text: 'Energy-loss spectrum' }), h('span', { class: 'spec-pos' })]), bCv]), readout])]),
    h('div', { class: 'ctrls' }, [sWin, sLog, sT, sR, sE, sView]));
  const posLab = root.querySelector('.spec-pos');
  function upd() { const n = MX * MY; status.textContent = mapping ? 'Mapping… ' + Math.round(acq / n * 100) + '%' : acq >= n ? MX + '×' + MY + ' spectrum image' : acq ? 'Paused' : ''; mapBtn.textContent = mapping ? 'Pause' : acq >= n ? 'Re-acquire map' : acq ? 'Resume map' : 'Acquire spectrum image'; }
  function changed() { acquire(); if (acq) { acq = 0; mapping = true; } paintMap(); upd(); }
  function acquire() {
    const c = eelsComp(probe[0], probe[1]), W = WIN[win], a = eelsAcquire(c, P, dose, W.e0, W.e1, W.d), fits = [];
    if (!W.log) for (const el of W.edges) { const r = edgeSignal(el, a.E, a.Y); if (r.A) fits.push({ el, A: r.A, r: r.r, sig: r.sig }); }
    spec = { E: a.E, Y: a.Y, fits, c };
    readout.innerHTML = '';
    if (W.log) readout.append(h('p', { text: 'Zero-loss peak, plasmon peaks at multiples of Eₚ (relative areas follow a Poisson distribution in t/λ), and the Li K edge at 55 eV.' }));
    else { readout.append(h('p', { class: 'key' }, [h('i', { class: 'k-fit' }), document.createTextNode('pre-edge fit window  '), h('i', { class: 'k-sig' }), document.createTextNode('background-subtracted signal')])); }
    if (!W.log) fits.forEach(f => readout.append(h('p', {}, [h('strong', { text: EDGE[f.el].lab + ' ' }), document.createTextNode('signal ' + sci(Math.max(0, f.sig)) + ' counts above a fitted A·E⁻ʳ background, r = ' + f.r.toFixed(2))])));
    posLab.textContent = 'x = ' + probe[0].toFixed(0) + ' nm';
    paintOverview(overview, c, P); dirty = true;
  }
  function paintMap() {
    const d = mImg.data, n = MX * MY, mx = {};
    MAP_EL.forEach(([el]) => { let m = 0; for (let i = 0; i < acq; i++) m = Math.max(m, maps[el][i]); mx[el] = m || 1; });
    for (let i = 0; i < n; i++) { if (i >= acq) { d[i * 4 + 3] = 0; continue; } let r = 0, g = 0, b = 0; for (const [el, c] of MAP_EL) { if (!show[el]) continue; const v = clamp(maps[el][i] / mx[el], 0, 1); r += c[0] * v * 0.8; g += c[1] * v * 0.8; b += c[2] * v * 0.8; } d[i * 4] = Math.min(255, r); d[i * 4 + 1] = Math.min(255, g); d[i * 4 + 2] = Math.min(255, b); d[i * 4 + 3] = 240; }
    mctx.putImageData(mImg, 0, 0); compose();
  }
  function compose() { cctx.clearRect(0, 0, 240, 120); cctx.drawImage(haadf, 0, 0); if (view === 'map' && acq) { cctx.imageSmoothingEnabled = true; cctx.drawImage(mapCv, 0, 0, 240, 120); } dirty = true; }
  function stepMap() {
    const t0 = performance.now(), n = MX * MY;
    while (acq < n && performance.now() - t0 < 10) {
      const x = (acq % MX + 0.5) * EX / MX, y = ((acq / MX | 0) + 0.5) * EY / MY, c = eelsComp(x, y);
      for (const [el, , d] of MAP_EL) { const w = FITW[el], a = eelsAcquire(c, P, dose, w[0], w[3], d); maps[el][acq] = Math.max(0, edgeSignal(el, a.E, a.Y).sig); }
      acq++;
    }
    if (acq >= n) mapping = false; paintMap(); upd();
  }
  function drawA() {
    const { ctx, w, h: hh } = PA; PA.reset(); ctx.imageSmoothingEnabled = true; ctx.drawImage(comp, 0, 0, w, hh);
    const sx = w / EX, X = probe[0] * sx, Y = probe[1] * sx;
    ctx.strokeStyle = COL.acc; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.arc(X, Y, 4, 0, TAU); ctx.moveTo(X - 10, Y); ctx.lineTo(X - 5, Y); ctx.moveTo(X + 5, Y); ctx.lineTo(X + 10, Y); ctx.moveTo(X, Y - 10); ctx.lineTo(X, Y - 5); ctx.moveTo(X, Y + 5); ctx.lineTo(X, Y + 10); ctx.stroke();
    scaleBar(ctx, w - 14 - 20 * sx, w - 14, hh - 12, '20 nm');
  }
  function drawB() {
    const { ctx, w, h: hh } = PB, W = WIN[win], L = 54, R = 12, T = 14, B = 32, pw = w - L - R, ph = hh - T - B; PB.reset(); ctx.fillStyle = COL.panel; ctx.fillRect(0, 0, w, hh);
    if (!spec) return; const { E, Y } = spec;
    let mx = 1, mnp = Infinity; for (let k = 0; k < Y.length; k++) { if (Y[k] > mx) mx = Y[k]; if (Y[k] > 0 && Y[k] < mnp) mnp = Y[k]; }
    const lmax = Math.log10(mx) + 0.15, lmin = Math.max(Math.floor(Math.log10(Math.max(mnp, 1))), lmax - 9);
    const xE = e => L + (e - W.e0) / (W.e1 - W.e0) * pw, yV = v => logY ? T + ph * (1 - (Math.log10(Math.max(v, Math.pow(10, lmin))) - lmin) / (lmax - lmin)) : T + ph * (1 - v / (mx * 1.06));
    ctx.strokeStyle = 'rgba(255,255,255,.08)'; ctx.lineWidth = 1; ctx.font = '500 10px ' + FONT; ctx.fillStyle = COL.t2;
    if (logY) { for (let e = Math.ceil(lmin); e <= lmax; e += Math.max(1, Math.round((lmax - lmin) / 5))) { const y = yV(Math.pow(10, e)); ctx.beginPath(); ctx.moveTo(L, y); ctx.lineTo(L + pw, y); ctx.stroke(); ctx.textAlign = 'right'; ctx.textBaseline = 'middle'; ctx.fillText('10' + String(e).split('').map(q => SUP[q]).join(''), L - 6, y); } }
    else { for (let i = 0; i <= 4; i++) { const v = mx * 1.06 * i / 4, y = yV(v); ctx.beginPath(); ctx.moveTo(L, y); ctx.lineTo(L + pw, y); ctx.stroke(); ctx.textAlign = 'right'; ctx.textBaseline = 'middle'; ctx.fillText(sci(v), L - 6, y); } }
    const span = W.e1 - W.e0, stepE = span > 300 ? 100 : span > 150 ? 50 : span > 60 ? 20 : 10;
    ctx.textAlign = 'center'; ctx.textBaseline = 'top'; for (let e = Math.ceil(W.e0 / stepE) * stepE; e <= W.e1; e += stepE) { const x = xE(e); ctx.fillText(String(e), x, T + ph + 6); ctx.beginPath(); ctx.moveTo(x, T + ph); ctx.lineTo(x, T + ph + 4); ctx.stroke(); }
    txt(ctx, 'energy loss (eV)', L + pw / 2, hh - 6, { align: 'center', size: 10 }); ctx.save(); ctx.translate(12, T + ph / 2); ctx.rotate(-Math.PI / 2); txt(ctx, 'counts', 0, 0, { align: 'center', size: 10 }); ctx.restore();
    ctx.save(); ctx.beginPath(); ctx.rect(L, T - 2, pw, ph + 2); ctx.clip();
    for (const f of spec.fits) {
      const w4 = FITW[f.el]; ctx.fillStyle = 'rgba(74,158,255,.10)'; ctx.fillRect(xE(w4[0]), T, xE(w4[1]) - xE(w4[0]), ph);
      ctx.beginPath(); let first = true; for (let k = 0; k < E.length; k++) { if (E[k] < w4[2] || E[k] > w4[3]) continue; const x = xE(E[k]), y = yV(Y[k]); if (first) { ctx.moveTo(x, y); first = false; } else ctx.lineTo(x, y); }
      for (let k = E.length - 1; k >= 0; k--) { if (E[k] < w4[2] || E[k] > w4[3]) continue; ctx.lineTo(xE(E[k]), yV(f.A * Math.pow(E[k], -f.r))); }
      ctx.closePath(); ctx.fillStyle = 'rgba(88,185,71,.32)'; ctx.fill();
    }
    ctx.beginPath(); for (let k = 0; k < E.length; k++) { const x = xE(E[k]), y = yV(Y[k]); if (k) ctx.lineTo(x, y); else ctx.moveTo(x, y); } ctx.strokeStyle = COL.t1; ctx.lineWidth = 1.1; ctx.stroke();
    for (const f of spec.fits) { const w4 = FITW[f.el]; ctx.setLineDash([5, 4]); ctx.strokeStyle = COL.blu; ctx.lineWidth = 1.5; ctx.beginPath(); for (let e = w4[0], i = 0; e <= w4[3]; e += (w4[3] - w4[0]) / 80, i++) { const x = xE(e), y = yV(f.A * Math.pow(e, -f.r)); if (i) ctx.lineTo(x, y); else ctx.moveTo(x, y); } ctx.stroke(); ctx.setLineDash([]); }
    ctx.restore();
    for (const el of W.edges) { const x = xE(EDGE[el].E0); ctx.strokeStyle = 'rgba(245,184,61,.6)'; ctx.setLineDash([2, 3]); ctx.beginPath(); ctx.moveTo(x, T); ctx.lineTo(x, T + ph); ctx.stroke(); ctx.setLineDash([]); txt(ctx, EDGE[el].lab, x + 4, T + 8, { color: COL.amb, size: 10, w: 600 }); }
    if (W.log) { const ep = spec.c.Ep; txt(ctx, 'ZLP', xE(0) + 6, T + 8, { color: COL.t1, size: 10, w: 600 }); txt(ctx, 'plasmon', xE(ep), T + 22, { color: COL.t1, size: 10, w: 600, align: 'center' }); }
  }
  function setProbe(x, y) { probe = [clamp(x, 0.5, EX - 0.5), clamp(y, 0.5, EY - 0.5)]; acquire(); }
  drag(PA, p => setProbe(p[0] / PA.w * EX, p[1] / PA.h * EY));
  drag(SC, (p, ph) => { const q = SC.cam.unproject(p[0], p[1], 0); if (ph === 'down' && (Math.abs(q[0]) > EX / 2 || Math.abs(q[1]) > EY / 2)) return; setProbe(q[0] + EX / 2, EY / 2 - q[1]); });
  aCv.addEventListener('keydown', e => { const m = { ArrowLeft: [-1, 0], ArrowRight: [1, 0], ArrowUp: [0, -1], ArrowDown: [0, 1] }[e.key]; if (m) { e.preventDefault(); setProbe(probe[0] + m[0], probe[1] + m[1]); } });
  whenVisible(root, () => {
    if (!started) { started = true; haadf = eelsTexture(); compose(); acquire(); upd(); }
    if (mapping) stepMap();
    if (dirty) { drawEelsScene(SC, { probe, tex: comp, specCv: overview, labels: true }); drawA(); drawB(); dirty = false; }
  });
}

/* ------------------------------------------------------------------ boot */
/* ======================================== LLZO 4D-STEM (precomputed with abTEM) */
/* Data: js/llzo-4dstem-data.js (window.LLZO_4DSTEM), multislice simulations of
   tetragonal Li7La3Zr2O12 (Materials Project mp-942733) along [001]. */
const LLZO_EL = { La: [245, 184, 61, 3.2], Zr: [63, 208, 224, 3.0], O: [255, 107, 90, 1.9], Li: [120, 220, 110, 1.7] };
function loadLLZO() {
  if (window.LLZO_4DSTEM) return Promise.resolve(window.LLZO_4DSTEM);
  return new Promise((res, rej) => { const s = document.createElement('script'); s.src = 'js/llzo-4dstem-data.js'; s.onload = () => (window.LLZO_4DSTEM ? res(window.LLZO_4DSTEM) : rej(new Error('no data'))); s.onerror = () => rej(new Error('missing data file')); document.head.appendChild(s); });
}
async function decodePreset(p) {
  if (p.I) return p.I;
  const bin = atob(p.b64), u = new Uint8Array(bin.length); for (let i = 0; i < bin.length; i++) u[i] = bin.charCodeAt(i);
  const buf = await new Response(new Blob([u]).stream().pipeThrough(new DecompressionStream('gzip'))).arrayBuffer();
  const v = new Uint8Array(buf), lut = new Float32Array(256), I = new Float32Array(v.length);
  for (let i = 1; i < 256; i++) lut[i] = p.M * Math.pow(10, p.decades * (i / 255 - 1));
  for (let i = 0; i < v.length; i++) I[i] = lut[v[i]];
  return (p.I = I);
}
function llzoWidget(root) {
  const sceneCv = h('canvas', { class: 'scene', role: 'img', 'aria-label': '3D view: electron probe over an LLZO crystal viewed along [001], with the simulated diffraction pattern on the detector. Drag to move the probe.' });
  const aCv = h('canvas', { tabindex: '0', 'aria-label': 'Virtual image of LLZO. Drag, or use the arrow keys, to move the probe.' }), bCv = h('canvas', { 'aria-label': 'Simulated diffraction pattern. Click a spot to make a dark-field image.' });
  const SC = makeCanvas(sceneCv, 760, 420), PA = makeCanvas(aCv, 360, 360), PB = makeCanvas(bCv, 360, 360);
  SC.cam = new Cam(-32, 28); stemCamFit(SC.cam, 760, 420, 16);
  const status = h('span', { class: 'hint', 'aria-live': 'polite', text: 'Loading simulated data…' }), aTitle = h('span'), contrast = h('span', { class: 'hint', style: 'color:var(--t1)' });
  let D = null, P = null, I = null, NS = 24, N = 64, mode = 'bf', atomsOn = true, gamma = 0.3, ix = 6, iy = 6, vdf = null, dirty = true, busy = false, vals = new Float32Array(NS * NS), mask = new Int32Array(0);
  const vCv = offscreen(2 * NS, 2 * NS), vctx = vCv.getContext('2d'), vImg = vctx.createImageData(2 * NS, 2 * NS);
  const dpCv = offscreen(N, N), dctx = dpCv.getContext('2d'), dImg = dctx.createImageData(N, N), TW = 384, comp = offscreen(TW, TW), cctx = comp.getContext('2d'), atomCv = offscreen(TW, TW);
  const TITLES = { bf: 'Virtual bright-field image', abf: 'Virtual annular bright-field image', adf: 'Virtual annular dark-field image', vdf: 'Virtual dark-field image' };
  const sPre = seg({ label: 'Convergence semi-angle', options: [['0', '1 mrad · nanobeam'], ['1', '6 mrad'], ['2', '20 mrad · atomic resolution']], value: '2', on: v => choose(+v) });
  const sDet = seg({ label: 'Virtual detector', options: [['bf', 'Bright field'], ['abf', 'Annular BF'], ['adf', 'Annular DF'], ['vdf', 'Dark field (spot)']], value: mode, on: v => { mode = v; if (v === 'vdf' && !vdf) vdf = brightest(); remask(); } });
  const sAt = seg({ label: 'Atom columns', options: [['on', 'Show'], ['off', 'Hide']], value: 'on', on: v => { atomsOn = v === 'on'; compose(); } });
  const sG = slider({ label: 'Display gamma', min: 0.1, max: 1, step: 0.05, value: gamma, fmt: v => v.toFixed(2), on: v => { gamma = v; paintDP(); } });
  const legend = h('div', { class: 'legend' }, Object.keys(LLZO_EL).map(k => h('span', { class: 'lg', style: '--c:rgb(' + LLZO_EL[k].slice(0, 3) + ')', text: k })));
  root.append(
    h('div', { class: 'sim-scene' }, sceneCv),
    h('div', { class: 'sim-row' }, [
      h('div', { class: 'panel' }, [h('div', { class: 'panel-h' }, [aTitle, h('span', { text: 'drag to move probe' })]), aCv, legend]),
      h('div', { class: 'panel' }, [h('div', { class: 'panel-h' }, [h('span', { text: 'Simulated diffraction pattern' }), h('span', { text: 'click a spot → dark field' })]), bCv])]),
    h('div', { class: 'sim-bar' }, [contrast, status]),
    h('div', { class: 'ctrls' }, [sPre, sDet, sAt, sG]));
  const A2 = () => 2 * D.a; // field of view: 2 x 2 unit cells
  const kOf = j => (j - P.c0) * P.dk;
  function titles() { aTitle.textContent = TITLES[mode] + ' · 2 × 2 unit cells'; }
  function brightest() {
    let best = -1, bx = 0, by = 0; const o = (iy * NS + ix) * N * N, r0 = P.alpha + 2 * P.dk;
    for (let y = 0; y < N; y++) for (let x = 0; x < N; x++) { const k = Math.hypot(kOf(x), kOf(y)); if (k < r0 || k > 0.9 * P.c0 * P.dk) continue; const v = I[o + y * N + x]; if (v > best) { best = v; bx = x; by = y; } }
    return centroid(bx, by);
  }
  function adfRange() { const kmax = (N - 1 - P.c0) * P.dk; return [Math.min(P.alpha * 1.1 + P.dk, kmax * 0.8), kmax]; }
  function remask() {
    const idx = [], ar = adfRange(), vr = Math.max(0.04, 2.2 * P.dk);
    for (let y = 0; y < N; y++) for (let x = 0; x < N; x++) {
      const kx = kOf(x), ky = kOf(y), k = Math.hypot(kx, ky);
      const m = mode === 'bf' ? k <= P.alpha : mode === 'abf' ? k >= P.alpha / 2 && k <= P.alpha : mode === 'adf' ? k >= ar[0] && k <= ar[1] : vdf && Math.hypot(kx - vdf[0], ky - vdf[1]) <= vr;
      if (m) idx.push(y * N + x);
    }
    mask = Int32Array.from(idx);
    let mn = Infinity, mx = -Infinity;
    for (let p = 0; p < NS * NS; p++) { let s = 0; const o = p * N * N; for (let j = 0; j < mask.length; j++) s += I[o + mask[j]]; vals[p] = s; if (s < mn) mn = s; if (s > mx) mx = s; }
    // never stretch less than a 10% span: a nearly uniform image should look uniform
    const span = Math.max(mx - mn, 0.1 * mx), lo = (mx + mn) / 2 - span / 2, sc = span > 0 ? 255 / span : 0, d = vImg.data;
    for (let y = 0; y < 2 * NS; y++) for (let x = 0; x < 2 * NS; x++) { const g = (vals[(y % NS) * NS + (x % NS)] - lo) * sc, q = (y * 2 * NS + x) * 4; d[q] = d[q + 1] = d[q + 2] = g; d[q + 3] = 255; }
    contrast.textContent = 'Image contrast (max − min) / max: ' + (mx > 0 ? (100 * (mx - mn) / mx).toFixed(mx - mn < 0.01 * mx ? 2 : 0) : '0') + '%';
    vctx.putImageData(vImg, 0, 0); titles(); compose();
  }
  function drawAtoms() {
    const c = atomCv.getContext('2d'), s = TW / A2(); c.clearRect(0, 0, TW, TW);
    for (const el of ['Li', 'O', 'Zr', 'La']) { const xy = D.atoms[el], e = LLZO_EL[el]; c.fillStyle = 'rgba(' + e.slice(0, 3) + ',.9)'; c.strokeStyle = 'rgba(0,0,0,.55)'; c.lineWidth = 0.8;
      for (let ty = 0; ty < 2; ty++) for (let tx = 0; tx < 2; tx++) for (let i = 0; i < xy.length; i += 2) { c.beginPath(); c.arc((xy[i] + tx * D.a) * s, (xy[i + 1] + ty * D.a) * s, e[3], 0, TAU); c.fill(); c.stroke(); } }
  }
  function compose() { cctx.fillStyle = '#0d1016'; cctx.fillRect(0, 0, TW, TW); cctx.imageSmoothingEnabled = true; cctx.drawImage(vCv, 0, 0, TW, TW); if (atomsOn) { cctx.globalAlpha = 0.85; cctx.drawImage(atomCv, 0, 0); cctx.globalAlpha = 1; } dirty = true; }
  function paintDP() {
    if (!I) return; const o = (iy * NS + ix) * N * N; let mx = 0; for (let i = 0; i < N * N; i++) if (I[o + i] > mx) mx = I[o + i];
    const d = dImg.data; for (let i = 0; i < N * N; i++) { const li = Math.min(255, (Math.pow(I[o + i] / mx, gamma) * 255) | 0) * 3; d[i * 4] = LUT[li]; d[i * 4 + 1] = LUT[li + 1]; d[i * 4 + 2] = LUT[li + 2]; d[i * 4 + 3] = 255; }
    dctx.putImageData(dImg, 0, 0); dirty = true;
  }
  async function choose(k) {
    if (busy) return; busy = true; status.textContent = 'Decoding…';
    try { P = D.presets[k]; I = await decodePreset(P); vdf = null; if (mode === 'vdf') vdf = brightest(); remask(); paintDP(); SC.map.kmax = P.c0 * P.dk;
      status.textContent = '· abTEM multislice, 300 kV, ' + (D.thickness / 10).toFixed(1) + ' nm thick, α = ' + P.mrad + ' mrad, ' + NS + '×' + NS + ' probe positions per unit cell';
    } catch (e) { console.error(e); status.textContent = 'Could not decode the data in this browser.'; }
    busy = false;
  }
  function setProbe(x, y) { // x, y in Å inside the 2x2-cell field
    const st = D.a / NS; ix = ((Math.round(x / st) % NS) + NS) % NS; iy = ((Math.round(y / st) % NS) + NS) % NS;
    probeA = [clamp(x, 0, A2()), clamp(y, 0, A2())]; paintDP();
  }
  let probeA = [0, 0];
  SC.map = { world: p => [p[0] / A2() * 2 * ST.H - ST.H, ST.H - p[1] / A2() * 2 * ST.H], kmax: 1, scanH: null };
  function drawA() {
    const { ctx, w } = PA; PA.reset(); ctx.drawImage(comp, 0, 0, w, w);
    const s = w / A2(), X = probeA[0] * s, Y = probeA[1] * s, rp = Math.max(3, 0.51 / P.alpha * s / 2);
    ctx.strokeStyle = COL.acc; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.arc(X, Y, rp, 0, TAU); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(X - rp - 6, Y); ctx.lineTo(X - rp - 1, Y); ctx.moveTo(X + rp + 1, Y); ctx.lineTo(X + rp + 6, Y); ctx.moveTo(X, Y - rp - 6); ctx.lineTo(X, Y - rp - 1); ctx.moveTo(X, Y + rp + 1); ctx.lineTo(X, Y + rp + 6); ctx.stroke();
    scaleBar(ctx, w - 14 - 5 * s, w - 14, w - 14, '5 Å');
  }
  function drawB() {
    const { ctx, w } = PB; PB.reset(); ctx.imageSmoothingEnabled = true; ctx.drawImage(dpCv, 0, 0, w, w);
    const s = w / (N * P.dk), cx = (P.c0 + 0.5) * w / N; ctx.strokeStyle = 'rgba(120,235,100,.95)'; ctx.lineWidth = 1.5;
    const ring = (x, y, r) => { ctx.beginPath(); ctx.arc(x, y, r, 0, TAU); ctx.stroke(); };
    if (mode === 'bf') ring(cx, cx, P.alpha * s);
    else if (mode === 'abf' || mode === 'adf') { const r = mode === 'abf' ? [P.alpha / 2, P.alpha] : adfRange(); ctx.fillStyle = 'rgba(88,185,71,.14)'; ctx.beginPath(); ctx.arc(cx, cx, r[1] * s, 0, TAU); ctx.arc(cx, cx, r[0] * s, 0, TAU, true); ctx.fill(); ring(cx, cx, r[0] * s); ring(cx, cx, r[1] * s); }
    else if (vdf) ring(cx + vdf[0] * s, cx + vdf[1] * s, Math.max(0.04, 2.2 * P.dk) * s);
    const bar = P.dk * N > 2 ? 1 : 0.5; scaleBar(ctx, w - 14 - bar * s, w - 14, w - 14, bar + ' Å⁻¹');
  }
  drag(PA, p => { if (P) setProbe(p[0] / PA.w * A2(), p[1] / PA.h * A2()); });
  drag(SC, (p, ph) => { if (!P) return; const q = SC.cam.unproject(p[0], p[1], 0); if (ph === 'down' && (Math.abs(q[0]) > ST.H || Math.abs(q[1]) > ST.H)) return; setProbe((q[0] + ST.H) / (2 * ST.H) * A2(), (ST.H - q[1]) / (2 * ST.H) * A2()); });
  aCv.addEventListener('keydown', e => { if (!P) return; const st = D.a / NS, m = { ArrowLeft: [-st, 0], ArrowRight: [st, 0], ArrowUp: [0, -st], ArrowDown: [0, st] }[e.key]; if (m) { e.preventDefault(); setProbe(probeA[0] + m[0], probeA[1] + m[1]); } });
  bCv.addEventListener('click', e => {
    if (!P) return; const p = pos(PB, e); let x = Math.floor(p[0] / PB.w * N), y = Math.floor(p[1] / PB.w * N), best = -1, bx = x, by = y; const o = (iy * NS + ix) * N * N;
    for (let yy = y - 3; yy <= y + 3; yy++) for (let xx = x - 3; xx <= x + 3; xx++) { if (xx < 0 || yy < 0 || xx >= N || yy >= N) continue; const v = I[o + yy * N + xx]; if (v > best) { best = v; bx = xx; by = yy; } }
    vdf = centroid(bx, by); mode = 'vdf'; sDet.set('vdf'); remask();
  });
  function centroid(bx, by) { // intensity-weighted centre of the disk around (bx, by)
    const o = (iy * NS + ix) * N * N, pk = I[o + by * N + bx]; let sx = 0, sy = 0, sw = 0;
    for (let yy = by - 4; yy <= by + 4; yy++) for (let xx = bx - 4; xx <= bx + 4; xx++) { if (xx < 0 || yy < 0 || xx >= N || yy >= N) continue; const v = I[o + yy * N + xx]; if (v < 0.3 * pk) continue; sx += v * xx; sy += v * yy; sw += v; }
    return sw ? [kOf(sx / sw), kOf(sy / sw)] : [kOf(bx), kOf(by)];
  }
  let started = false;
  whenVisible(root, () => {
    if (!started) {
      started = true;
      loadLLZO().then(d => { D = d; NS = d.ns; N = d.n; drawAtoms(); setProbe(D.a * 0.5, D.a * 0.5); return choose(2); })
        .catch(() => { status.textContent = 'The simulation data file (js/llzo-4dstem-data.js) was not found. Keep the js folder next to index.html.'; });
    }
    if (dirty && P && I) { drawStemScene(SC, { alpha: P.alpha, df: 0 }, comp, dpCv, probeA, true); drawA(); drawB(); dirty = false; }
  });
}

function boot() {
  const P = { fib: fibPillar, stem: stemPillar, eels: eelsPillar }, Wd = { fib: fibWidget, stem: stemWidget, eels: eelsWidget, llzo: llzoWidget };
  document.querySelectorAll('canvas[data-pillar]').forEach(cv => { try { P[cv.dataset.pillar](cv); } catch (e) { console.error(e); } });
  document.querySelectorAll('[data-sim]').forEach(root => { try { Wd[root.dataset.sim](root); } catch (e) { fail(root, e); } });
}
if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', boot); else boot();
})();

