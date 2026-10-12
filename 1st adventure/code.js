// Quantum Plants vs. Zombies: Adventure 1, "Qubit basics" (3 levels)
// A. One qubit: X, H, Z, measurement in the Z or X basis, classical vs quantum zombies.
// B. Adds RX, RY, RZ. Quantum zombies carry visible superposed states; a qubit that meets one is
//    tested against that state by the Born rule, so the rotation angles must reproduce it.
// C. Two qubits + CNOT. The boss is a Bell-state quantum zombie. A product state overlaps a Bell
//    state with probability at most 1/2, so only an entangled register can wear it down.
// Bit strings read q1 q0, as in Qiskit.
(() => {
  "use strict";

  // ==QUANTUM CORE==
  let NQ = 1, DIM = 2;
  const setQubits = n => { NQ = n; DIM = 1 << n; };
  const SQ = Math.SQRT1_2;
  const C = (re, im = 0) => [re, im];
  const cadd = (a, b) => [a[0] + b[0], a[1] + b[1]];
  const csub = (a, b) => [a[0] - b[0], a[1] - b[1]];
  const cmul = (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]];
  const cconj = a => [a[0], -a[1]];
  const cabs2 = a => a[0] * a[0] + a[1] * a[1];
  const bitOf = (i, q) => (i >> q) & 1;
  const zeros = n => Array.from({ length: n }, () => C(0));
  const basisState = i => { const a = zeros(DIM); a[i] = C(1); return a; };
  function normalize(a) {
    const n = Math.sqrt(a.reduce((s, v) => s + cabs2(v), 0)) || 1;
    return a.map(v => [v[0] / n, v[1] / n]);
  }
  const inner = (u, v) => u.reduce((s, x, i) => cadd(s, cmul(cconj(x), v[i])), C(0));

  const MAT = {
    X: () => [[C(0), C(1)], [C(1), C(0)]],
    H: () => [[C(SQ), C(SQ)], [C(SQ), C(-SQ)]],
    Z: () => [[C(1), C(0)], [C(0), C(-1)]],
    RX: t => { const c = Math.cos(t / 2), s = Math.sin(t / 2); return [[C(c), C(0, -s)], [C(0, -s), C(c)]]; },
    RY: t => { const c = Math.cos(t / 2), s = Math.sin(t / 2); return [[C(c), C(-s)], [C(s), C(c)]]; },
    RZ: t => { const c = Math.cos(t / 2), s = Math.sin(t / 2); return [[C(c, -s), C(0)], [C(0), C(c, s)]]; }
  };
  function apply1(a, q, m) {
    const b = a.map(v => v.slice());
    for (let i = 0; i < a.length; i++) {
      if (bitOf(i, q)) continue;
      const j = i | (1 << q);
      b[i] = cadd(cmul(m[0][0], a[i]), cmul(m[0][1], a[j]));
      b[j] = cadd(cmul(m[1][0], a[i]), cmul(m[1][1], a[j]));
    }
    return b;
  }
  const applyCNOT = (a, c, t) => a.map((v, i) => (bitOf(i, c) ? a[i ^ (1 << t)] : v).slice());
  // A circuit column holds one single-qubit gate per lane, or one CNOT.
  function applyCell(a, cell) {
    if (!cell) return a;
    if (cell.cx) return applyCNOT(a, cell.cx.c, cell.cx.t);
    let b = a;
    cell.g.forEach((g, q) => { if (g) b = apply1(b, q, MAT[g.type](g.theta)); });
    return b;
  }
  const runCells = cells => cells.reduce((a, c) => applyCell(a, c), basisState(0));
  const probs = a => a.map(cabs2);
  function sampleIndex(p) {
    let r = Math.random() * p.reduce((s, v) => s + v, 0);
    for (let i = 0; i < p.length; i++) { r -= p[i]; if (r < 0) return i; }
    return p.length - 1;
  }
  // Measure every qubit in the Z or X basis (one joint sample), returning the collapsed state.
  function measureAll(a, basis) {
    const toX = s => { let b = s; for (let q = 0; q < NQ; q++) b = apply1(b, q, MAT.H()); return b; };
    const k = sampleIndex(probs(basis === "X" ? toX(a) : a));
    return { k, amps: basis === "X" ? toX(basisState(k)) : basisState(k) };
  }
  // Component of the register along |psi⟩ on qubit q: pairs (i, j, ⟨psi|_q a⟩).
  function qubitComponent(a, q, psi) {
    const comp = [];
    let p = 0;
    for (let i = 0; i < a.length; i++) {
      if (bitOf(i, q)) continue;
      const j = i | (1 << q);
      const v = cadd(cmul(cconj(psi[0]), a[i]), cmul(cconj(psi[1]), a[j]));
      comp.push([i, j, v]); p += cabs2(v);
    }
    return { comp, p };
  }
  const overlapQubit = (a, q, psi) => qubitComponent(a, q, psi).p;
  // Projective test of qubit q against |psi⟩: outcomes {|psi⟩⟨psi|, 1 − |psi⟩⟨psi|}, Born rule.
  function projectQubit(a, q, psi) {
    const { comp, p } = qubitComponent(a, q, psi);
    const hit = Math.random() < p;
    let b = hit ? zeros(a.length) : a.map(v => v.slice());
    for (const [i, j, v] of comp) {
      if (hit) { b[i] = cmul(psi[0], v); b[j] = cmul(psi[1], v); }
      else { b[i] = csub(b[i], cmul(psi[0], v)); b[j] = csub(b[j], cmul(psi[1], v)); }
    }
    return { hit, p, amps: normalize(b) };
  }
  // Projective test of the whole register against |phi⟩.
  function projectState(a, phi) {
    const v = inner(phi, a), p = cabs2(v);
    const hit = Math.random() < p;
    const b = hit ? phi.map(x => x.slice()) : a.map((x, i) => csub(x, cmul(phi[i], v)));
    return { hit, p, amps: normalize(b) };
  }
  // Bloch vector of qubit q (reduced state, so its length drops below 1 when entangled).
  function blochOf(a, q) {
    let z = 0, rho01 = C(0);
    for (let i = 0; i < a.length; i++) {
      z += bitOf(i, q) ? -cabs2(a[i]) : cabs2(a[i]);
      if (!bitOf(i, q)) rho01 = cadd(rho01, cmul(a[i], cconj(a[i | (1 << q)])));
    }
    const x = 2 * rho01[0], y = -2 * rho01[1];
    return { x, y, z, r: Math.hypot(x, y, z) };
  }
  const concurrence = a => 2 * Math.sqrt(cabs2(csub(cmul(a[0], a[3]), cmul(a[1], a[2]))));
  // cos(θ/2)|0⟩ + e^{iφ} sin(θ/2)|1⟩
  const qubitState = (th, ph) => [C(Math.cos(th / 2)), C(Math.sin(th / 2) * Math.cos(ph), Math.sin(th / 2) * Math.sin(ph))];
  const BELL = {
    "Φ+": [C(SQ), C(0), C(0), C(SQ)], "Φ−": [C(SQ), C(0), C(0), C(-SQ)],
    "Ψ+": [C(0), C(SQ), C(SQ), C(0)], "Ψ−": [C(0), C(SQ), C(-SQ), C(0)]
  };
  // ==END CORE==

  // ---------- State formatting ----------
  const PI = Math.PI;
  const gcd = (a, b) => (b ? gcd(b, a % b) : Math.abs(a));
  function piFrac(t) {
    const k = Math.round(t / (PI / 12));
    if (Math.abs(t - (k * PI) / 12) > 1e-6) return `${((t * 180) / PI).toFixed(1)}°`;
    if (k === 0) return "0";
    const g = gcd(k, 12), p = k / g, q = 12 / g;
    const num = `${p < 0 ? "−" : ""}${Math.abs(p) === 1 ? "" : Math.abs(p)}π`;
    return q === 1 ? num : `${num}/${q}`;
  }
  function piJS(t) {
    const k = Math.round(t / (PI / 12));
    if (Math.abs(t - (k * PI) / 12) > 1e-6) return t.toFixed(6);
    if (k === 0) return "0";
    const g = gcd(k, 12), p = k / g, q = 12 / g;
    const num = p === 1 ? "Math.PI" : p === -1 ? "-Math.PI" : `${p} * Math.PI`;
    return q === 1 ? num : `${num} / ${q}`;
  }
  const MAGS = [[1, "1"], [Math.sqrt(3) / 2, "√3/2"], [SQ, "1/√2"], [0.5, "1/2"]];
  const magStr = m => { const f = MAGS.find(([v]) => Math.abs(v - m) < 1e-6); return f ? f[1] : m.toFixed(3); };
  const wrapPhase = x => Math.atan2(Math.sin(x), Math.cos(x));
  const ketLabel = (i, n) => `|${i.toString(2).padStart(n, "0")}⟩`;

  // Writes a state as a ket sum, up to a global phase, e.g. "√3/2|0⟩ + i/2|1⟩" or "(|00⟩ + |11⟩)/√2".
  function fmtKet(a) {
    const n = Math.round(Math.log2(a.length));
    const k0 = a.findIndex(v => cabs2(v) > 1e-10);
    if (k0 < 0) return "0";
    const ph0 = Math.atan2(a[k0][1], a[k0][0]);
    const terms = [];
    a.forEach((v, i) => {
      const m = Math.hypot(v[0], v[1]);
      if (m > 1e-6) terms.push({ i, m, ph: wrapPhase(Math.atan2(v[1], v[0]) - ph0) });
    });
    const m0 = terms[0].m;
    const common = terms.length > 1 && terms.every(t => Math.abs(t.m - m0) < 1e-6) &&
      ((terms.length === 2 && Math.abs(m0 - SQ) < 1e-6) || (terms.length === 4 && Math.abs(m0 - 0.5) < 1e-6));
    const near = (ph, x) => Math.abs(wrapPhase(ph - x)) < 1e-6;
    const parts = terms.map((t, idx) => {
      let sign = "+", im = false, extra = "";
      if (near(t.ph, 0)) { /* real positive */ }
      else if (near(t.ph, PI)) sign = "−";
      else if (near(t.ph, PI / 2)) im = true;
      else if (near(t.ph, -PI / 2)) { sign = "−"; im = true; }
      else extra = `e^(i${piFrac(t.ph)})`;
      const mag = common ? "1" : magStr(t.m);
      let body;
      if (im) body = mag === "1" ? "i" : mag.startsWith("1/") ? `i/${mag.slice(2)}` : `i${mag}`;
      else if (extra) body = mag === "1" ? extra : `${mag}·${extra}`;
      else body = mag === "1" ? "" : mag;
      const ket = ketLabel(t.i, n);
      if (idx === 0) return `${sign === "−" ? "−" : ""}${body}${ket}`;
      return ` ${sign} ${body}${ket}`;
    }).join("");
    return common ? `(${parts})/${terms.length === 2 ? "√2" : "2"}` : parts;
  }
  const NAMED1 = [["|+⟩", qubitState(PI / 2, 0)], ["|−⟩", qubitState(PI / 2, PI)],
    ["|+i⟩", qubitState(PI / 2, PI / 2)], ["|−i⟩", qubitState(PI / 2, -PI / 2)],
    ["|0⟩", qubitState(0, 0)], ["|1⟩", qubitState(PI, 0)]];
  function nameOf(a) {
    const list = a.length === 2 ? NAMED1 : a.length === 4 ? Object.entries(BELL) : [];
    const f = list.find(([, v]) => cabs2(inner(v, a)) > 0.9999);
    return f ? f[0] : null;
  }
  function fmtNamed(a) {
    const k = fmtKet(a), nm = nameOf(a);
    return nm && nm !== k ? `${nm} = ${k}` : k;
  }
  // Short label for a qubit inside a register (drawn above the orb).
  function shortLabel(a, q) {
    const b = blochOf(a, q);
    if (b.r < 0.98) return "ent";
    const v = [[b.z, "|0⟩"], [-b.z, "|1⟩"], [b.x, "|+⟩"], [-b.x, "|−⟩"], [b.y, "|+i⟩"], [-b.y, "|−i⟩"]].find(([c]) => c > 0.98);
    return v ? v[1] : "ψ";
  }

  // ---------- Levels ----------
  const ROT = ["RX", "RY", "RZ"];
  const LEVELS = [
    { nq: 1, key: "A", name: "First qubit", slots: 4, tools: ["X", "H", "Z"],
      intro: "One qubit, one lane. Classical zombies (0, 1, +, −) fall to two matching measured outcomes beyond the gate. Quantum zombies walk through the gate and are only hurt by unmeasured qubits in the lab zone.",
      wave: [{ t: 0, type: "0" }, { t: 6, type: "0" }, { t: 13, type: "Q", st: "+" }, { t: 21, type: "1" },
        { t: 29, type: "+" }, { t: 37, type: "Q", st: "-" }, { t: 45, type: "-" }, { t: 52, type: "1" }, { t: 60, type: "Q", st: "+" }] },
    { nq: 1, key: "B", name: "Rotations", slots: 4, tools: ["X", "H", "Z", "RX", "RY", "RZ"],
      intro: "Quantum zombies now carry arbitrary superpositions, written on their cards. Set θ and plant RX, RY or RZ so the lane prepares exactly that state: the closer the match, the likelier each hit.",
      wave: [{ t: 0, type: "Q", st: "rand" }, { t: 9, type: "1" }, { t: 15, type: "Q", st: "rand" }, { t: 27, type: "Q", st: "rand" },
        { t: 35, type: "+" }, { t: 41, type: "Q", st: "rand" }, { t: 53, type: "Q", st: "rand" }, { t: 65, type: "Q", st: "rand" }] },
    { nq: 2, key: "C", name: "Entangled boss", slots: 5, tools: ["X", "H", "Z", "RX", "RY", "RZ", "CNOT"], boss: true,
      intro: "Two lanes, and a Bell-state boss. It passes the measurement gate and is tested against whole unmeasured registers. A product state matches a Bell state at most half the time, so you need entanglement to win.",
      wave: [{ t: 4, type: "0", lane: 0 }, { t: 12, type: "Q", st: "+", lane: 1 }, { t: 22, type: "1", lane: 1 },
        { t: 32, type: "-", lane: 0 }, { t: 44, type: "Q", st: "-", lane: 0 }] }
  ];
  const RAND_STATES = [];
  for (const th of [PI / 3, PI / 2, (2 * PI) / 3]) for (const ph of [0, PI / 2, PI, -PI / 2]) {
    if (th === PI / 2 && (ph === 0 || ph === PI)) continue;
    RAND_STATES.push([th, ph]);
  }
  function stateFor(st) {
    if (st === "+") return qubitState(PI / 2, 0);
    if (st === "-") return qubitState(PI / 2, PI);
    const [th, ph] = RAND_STATES[Math.floor(Math.random() * RAND_STATES.length)];
    return qubitState(th, ph);
  }

  const $ = id => document.getElementById(id);
  const canvas = $("lawn");
  const ctx = canvas.getContext("2d");
  const bcv = $("blochCv");
  const bctx = bcv.getContext("2d");

  const W = 1200, TOP = 44, FOOT = 52;
  const X0 = 100, XM = 584, XE = 1172;
  const FONT = '"Lilita One", "Arial Rounded MT Bold", "Trebuchet MS", system-ui, sans-serif';
  const LANE_END = 15, MEASURE_X = 5, QSPEED = 1, EMIT = 2, ZSPEED = 0.45, BOSS_SPEED = 0.2, BOSS_LAB_SPEED = 0.1, RAD = 0.3;
  const Q_HP = 2;
  // Lower end of the 3σ (Wilson) confidence interval for a success rate, from k successes in n tests.
  function wilsonLow(k, n, z = 3) {
    if (!n) return 0;
    const p = k / n, z2 = z * z;
    return (p + z2 / (2 * n) - z * Math.sqrt(p * (1 - p) / n + z2 / (4 * n * n))) / (1 + z2 / n);
  }
  let LANE_H = 200, DATA_H = 200, H = 296;
  let SLOTS = [];

  const xPos = u => (u <= MEASURE_X ? X0 + (u * (XM - X0)) / MEASURE_X
    : XM + ((u - MEASURE_X) * (XE - XM)) / (LANE_END - MEASURE_X));
  const laneY = q => TOP + q * LANE_H + LANE_H / 2;
  const lanes = () => Array.from({ length: NQ }, (_, q) => q);
  const emptyCell = () => ({ g: new Array(NQ).fill(null), cx: null });

  const G = {
    level: 0, cells: [], tool: "X", theta: PI / 3, basis: "Z", showFid: false,
    phase: "plan", baseHP: 10, gateHP: 10, score: 0,
    regs: [], zombies: [], boss: null, splats: [], wave: [], waveIdx: 0, waveT: 0,
    emitT: 0, paused: false, nextId: 1, gateFlash: 0, houseFlash: 0, monT: 0, lastOutcome: "", testLane: 0
  };
  const lvl = () => LEVELS[G.level];
  const gateAlive = () => G.gateHP > 0;
  const sym = s => (s === "-" ? "−" : s);

  function fitCanvas() {
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    canvas.width = W * dpr; canvas.height = H * dpr;
    canvas.style.aspectRatio = `${W} / ${H}`;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    bcv.width = 400 * dpr; bcv.height = 200 * dpr;
    bctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  }
  function say(t) { $("msg").textContent = t; }

  function loadLevel(i, keepStats) {
    G.level = i;
    const L = lvl();
    setQubits(L.nq);
    LANE_H = L.nq === 1 ? 200 : 170;
    DATA_H = LANE_H * L.nq;
    H = TOP + DATA_H + FOOT;
    SLOTS = Array.from({ length: L.slots }, (_, k) => ({ i: k, u: 0.85 + 0.8 * k }));
    G.cells = SLOTS.map(emptyCell);
    if (!keepStats) { G.baseHP = 10; G.score = 0; }
    G.gateHP = 10;
    G.basis = "Z";
    G.paused = false;
    G.splats = [];
    G.regs = [];
    G.zombies = [];
    G.waveIdx = 0; G.waveT = 0; G.emitT = 0;
    G.testLane = 0;
    G.wave = L.wave.map(w => {
      const psi = w.type === "Q" ? stateFor(w.st) : null;
      return { ...w, lane: w.lane ?? 0, psi, label: psi ? fmtNamed(psi) : sym(w.type), done: false };
    });
    G.boss = null;
    if (L.boss) {
      const keys = Object.keys(BELL), key = keys[Math.floor(Math.random() * keys.length)];
      G.boss = { key, phi: BELL[key], label: `${key} = ${fmtKet(BELL[key])}`, x: 13.6, prevX: 13.6,
        matches: 0, tests: 0, flash: 0, hitGate: false };
    }
    if (!L.tools.includes(G.tool)) G.tool = "X";
    G.phase = "plan";
    fitCanvas();
    buildTools();
    buildBars();
    syncControls();
    renderWave();
    say(`Level ${L.key}: ${L.name}. ${L.intro} Press Ready when your lane is planted.`);
    circuitChanged();
    autoTutorial();
  }

  function startWave() {
    if (G.phase !== "plan") return;
    G.phase = "battle";
    G.emitT = EMIT;
    say(`Level ${lvl().key} is moving. A fresh |${"0".repeat(NQ)}⟩ register leaves the base every 2 seconds.`);
    updateHUD();
  }

  function checkCleared() {
    if (G.phase !== "battle" || G.waveIdx < G.wave.length || G.zombies.length || G.boss) return;
    G.score += 25;
    G.phase = "cleared";
    const from = G.level, last = from === LEVELS.length - 1;
    say(last ? "Adventure 1 complete: superposition, measurement, phase and entanglement mastered. Level C restarts with a new boss."
      : `Level ${lvl().key} cleared. Level ${LEVELS[G.level + 1].key} is lining up.`);
    setTimeout(() => {
      if (G.phase !== "cleared" || G.level !== from) return;
      loadLevel(last ? G.level : G.level + 1, true);
    }, 1800);
  }

  // ---------- Entities ----------
  function emit() {
    G.regs.push({ id: G.nextId++, x: 0.25, prevX: 0.25, amps: basisState(0), next: 0, measured: false,
      sym: null, basis: null, lanes: new Array(NQ).fill(true), broken: false, flash: 0, met: new Set() });
  }
  function passSlots(r) {
    while (r.next < SLOTS.length && r.x >= SLOTS[r.next].u) {
      r.amps = applyCell(r.amps, G.cells[r.next]);
      r.next++;
    }
  }
  function measureReg(r) {
    const m = measureAll(r.amps, G.basis);
    r.amps = m.amps; r.measured = true; r.basis = G.basis; r.flash = 0.3;
    r.sym = lanes().map(q => (G.basis === "Z" ? String(bitOf(m.k, q)) : (bitOf(m.k, q) ? "-" : "+")));
    G.gateFlash = 0.25;
    G.lastOutcome = `${G.basis} basis → ${lanes().reverse().map(q => sym(r.sym[q])).join("")}`;
  }
  function spawnZombie(w) {
    const q = w.type === "Q";
    G.zombies.push({ id: G.nextId++, type: w.type, lane: Math.min(w.lane, NQ - 1), x: LANE_END, prevX: LANE_END,
      hits: 0, hp: q ? Q_HP : 0, psi: w.psi, label: w.label, alive: true, flash: 0, phase: Math.random() * 6 });
  }

  const crossed = (a0, a1, b0, b1, rad) => {
    const s = a0 - b0, e = a1 - b1;
    return Math.abs(s) <= rad || Math.abs(e) <= rad || s * e < 0;
  };
  function splat(u, q, color) { G.splats.push({ x: xPos(u), y: laneY(q) - 8, r: 6, life: 0.4, color }); }
  function defeat(z, pts, t) { z.alive = false; G.score += pts; splat(z.x, z.lane, "#b6f27a"); say(t); }
  function damageGate(n, t) {
    if (!gateAlive()) return;
    G.gateHP = Math.max(0, G.gateHP - n);
    G.gateFlash = 0.4;
    say(G.gateHP === 0
      ? "Measurement gate destroyed. Registers now stay quantum all the way to x = 15, and classical zombies walk to the base."
      : `${t}. Gate HP ${G.gateHP}.`);
  }
  function damageBase(n, t) {
    G.baseHP = Math.max(0, G.baseHP - n);
    G.houseFlash = 0.8;
    say(`${t}. Base HP ${G.baseHP}.`);
    if (G.baseHP === 0) { G.phase = "over"; say("Base lost. Restart to play again."); }
  }

  // ---------- Update ----------
  function update(dt) {
    G.gateFlash = Math.max(0, G.gateFlash - dt);
    G.houseFlash = Math.max(0, G.houseFlash - dt);
    for (const s of G.splats) { s.life -= dt; s.r += dt * 50; }
    G.splats = G.splats.filter(s => s.life > 0);
    if (G.phase !== "battle" || G.paused) return;

    G.waveT += dt;
    while (G.waveIdx < G.wave.length && G.waveT >= G.wave[G.waveIdx].t) {
      const w = G.wave[G.waveIdx++];
      w.done = true;
      spawnZombie(w);
      renderWave();
    }
    G.emitT += dt;
    while (G.emitT >= EMIT) { G.emitT -= EMIT; emit(); }

    for (const r of G.regs) {
      r.prevX = r.x; r.x += QSPEED * dt; r.flash = Math.max(0, r.flash - dt);
      passSlots(r);
      if (!r.measured && gateAlive() && r.x >= MEASURE_X) measureReg(r);
      if (r.x > LANE_END + 0.5 || r.lanes.every(l => !l)) r.dead = true;
    }

    for (const z of G.zombies) {
      if (!z.alive) continue;
      z.prevX = z.x; z.x -= ZSPEED * dt; z.flash = Math.max(0, z.flash - dt);
      if (z.type !== "Q" && gateAlive() && z.prevX >= MEASURE_X && z.x < MEASURE_X) {
        z.alive = false; splat(MEASURE_X, z.lane, "#ffd23f");
        damageGate(1, "A classical zombie smashed into the measurement gate");
        continue;
      }
      if (z.x <= 0.3) { z.alive = false; damageBase(1, `A ${z.type === "Q" ? "quantum" : "classical"} zombie reached the base`); }
    }

    const B = G.boss;
    if (B) {
      B.prevX = B.x; B.x -= (B.x < MEASURE_X ? BOSS_LAB_SPEED : BOSS_SPEED) * dt; B.flash = Math.max(0, B.flash - dt);
      if (!B.hitGate && B.x <= MEASURE_X) {
        B.hitGate = true;
        if (gateAlive()) damageGate(2, "The entangled boss shoved through the measurement gate. It is in the lab zone now");
      }
      if (B.x <= 0.3) {
        damageBase(3, "The boss reached the base and walks back to the far end (your matches so far still count)");
        Object.assign(B, { x: 13.6, prevX: 13.6, hitGate: false });
      }
    }

    // Measured registers against classical zombies (beyond the gate).
    for (const r of G.regs) {
      if (!r.measured) continue;
      for (const q of lanes()) {
        if (!r.lanes[q]) continue;
        for (const z of G.zombies) {
          if (!z.alive || z.type === "Q" || z.lane !== q) continue;
          const cx = (r.prevX + r.x + z.prevX + z.x) / 4;
          if (cx < MEASURE_X - 0.05 || !crossed(r.prevX, r.x, z.prevX, z.x, RAD)) continue;
          r.lanes[q] = false; r.broken = true;
          if (r.sym[q] === z.type) {
            z.hits++; z.flash = 0.2; splat(z.x, q, "#b6f27a");
            if (z.hits >= 2) defeat(z, 10, `Two matching ${sym(z.type)} outcomes defeated a classical zombie${NQ > 1 ? ` on q${q}` : ""}.`);
            else say(`Matching ${sym(z.type)}. One more hit needed.`);
          } else {
            splat(z.x, q, "#d9d9d9");
            say(`Outcome ${sym(r.sym[q])} does not match zombie ${sym(z.type)}.`);
          }
          break;
        }
      }
    }

    // Unmeasured registers against quantum zombies and the boss (lab zone, or the whole lane once the gate is gone).
    const end = gateAlive() ? MEASURE_X : LANE_END;
    // The boss only tests finished registers: ones that have passed every planted gate.
    const lastPlanted = G.cells.reduce((m, c, k) => (c.cx || c.g.some(Boolean) ? k : m), -1);
    for (const r of G.regs) {
      if (r.measured || r.dead) continue;
      if (G.boss && !r.broken && r.lanes.every(Boolean) && r.next > lastPlanted) {
        const Bs = G.boss;
        const cx = (r.prevX + r.x + Bs.prevX + Bs.x) / 4;
        if (cx >= SLOTS[0].u && cx <= end && !r.met.has("boss") && crossed(r.prevX, r.x, Bs.prevX, Bs.x, RAD)) {
          // One Bell-state test per register. A match uses the register up; a mismatch passes through.
          // The boss never heals: it falls once the matches prove a match rate above 1/2, the most any
          // un-entangled pair can reach (an entanglement witness).
          r.met.add("boss");
          const m = projectState(r.amps, Bs.phi);
          r.amps = m.amps;
          Bs.tests++;
          if (m.hit) {
            Bs.matches++; Bs.flash = 0.2;
            r.lanes = r.lanes.map(() => false); r.dead = true;
            lanes().forEach(q => splat(Bs.x, q, "#b6f27a"));
            say(`The register matched ${Bs.key} (chance ${(m.p * 100).toFixed(0)}%). Matches ${Bs.matches} of ${Bs.tests}.`);
          } else {
            lanes().forEach(q => splat(Bs.x, q, "rgba(201,167,255,.5)"));
            say(`The register did not match ${Bs.key} (chance was ${(m.p * 100).toFixed(0)}%) and passed through. Matches ${Bs.matches} of ${Bs.tests}.`);
          }
          if (wilsonLow(Bs.matches, Bs.tests) > 0.5) {
            G.score += 50;
            lanes().forEach(q => splat(Bs.x, q, "#ffd23f"));
            say(`Boss defeated. ${Bs.matches} matches in ${Bs.tests} tests prove a match rate above 50%, which no un-entangled pair can reach.`);
            G.boss = null;
          }
          if (r.dead) continue;
        }
      }
      for (const q of lanes()) {
        if (!r.lanes[q]) continue;
        for (const z of G.zombies) {
          if (!z.alive || z.type !== "Q" || z.lane !== q) continue;
          const cx = (r.prevX + r.x + z.prevX + z.x) / 4;
          if (cx < SLOTS[0].u || cx > end || r.met.has(z.id) || !crossed(r.prevX, r.x, z.prevX, z.x, RAD)) continue;
          // One projective test per qubit and zombie. A match uses the qubit up; a mismatch means no
          // interaction: the qubit passes through in its post-test state (unchanged if it was orthogonal).
          r.met.add(z.id);
          const m = projectQubit(r.amps, q, z.psi);
          r.amps = m.amps;
          const where = NQ > 1 ? `q${q} ` : "The qubit ";
          if (m.hit) {
            r.lanes[q] = false; r.broken = true;
            z.hp -= 1; z.flash = 0.2; splat(z.x, q, "#b6f27a");
            if (z.hp <= 0) defeat(z, 20, `${where}matched ${z.label} (chance ${(m.p * 100).toFixed(0)}%). Quantum zombie dissolved.`);
            else say(`${where}matched ${z.label} (chance ${(m.p * 100).toFixed(0)}%). One more match needed.`);
          } else {
            splat(z.x, q, "rgba(201,167,255,.5)");
            say(`${where}did not match ${z.label} (chance was ${(m.p * 100).toFixed(0)}%), so it passed straight through without interacting.`);
          }
          break;
        }
      }
    }

    G.regs = G.regs.filter(r => !r.dead);
    G.zombies = G.zombies.filter(z => z.alive);
    checkCleared();
  }

  // ---------- Circuit editing ----------
  function editCell(cell, lane) {
    const t = G.tool;
    if (t === "ERASE") { if (cell.cx) cell.cx = null; else cell.g[lane] = null; return true; }
    if (t === "CNOT") {
      const same = cell.cx && cell.cx.c === lane;
      cell.g = new Array(NQ).fill(null);
      cell.cx = same ? null : { c: lane, t: 1 - lane };
      return true;
    }
    if (cell.cx) cell.cx = null;
    const cur = cell.g[lane], rot = ROT.includes(t);
    const same = cur && cur.type === t && (!rot || Math.abs(cur.theta - G.theta) < 1e-9);
    cell.g[lane] = same ? null : { type: t, theta: rot ? G.theta : 0 };
    return true;
  }
  const gateText = g => (ROT.includes(g.type) ? `${g.type}(${piFrac(g.theta)})` : g.type);

  function qiskitCode() {
    const L = ["import { QuantumCircuit } from 'qiskit';", "",
      `const qc = QuantumCircuit(${NQ}, ${NQ});   // lane q = qubit q, all start in |0⟩`];
    G.cells.forEach((c, k) => {
      if (c.cx) { L.push(`qc.cx(${c.cx.c}, ${c.cx.t});   // slot ${k + 1}`); return; }
      c.g.forEach((g, q) => {
        if (!g) return;
        const n = g.type.toLowerCase();
        L.push(ROT.includes(g.type) ? `qc.${n}(${piJS(g.theta)}, ${q});   // slot ${k + 1}` : `qc.${n}(${q});   // slot ${k + 1}`);
      });
    });
    const all = lanes().join(", ");
    if (G.basis === "X") L.push(`${lanes().map(q => `qc.h(${q});`).join(" ")}   // X-basis readout`);
    L.push(`qc.measure([${all}], [${all}]);`,
      "// next step: run on the pre-created Qollab backend; keys read high qubit first.");
    return L.join("\n");
  }

  function currentTarget() {
    if (G.boss) return { kind: "boss", vec: G.boss.phi, label: `Boss ${G.boss.label}` };
    const live = G.zombies.filter(z => z.type === "Q").sort((a, b) => a.x - b.x)[0];
    if (live) return { kind: "q", vec: live.psi, lane: live.lane, label: `${live.label}${NQ > 1 ? ` on q${live.lane}` : ""}` };
    const next = G.wave.find(w => !w.done && w.type === "Q");
    if (next) return { kind: "q", vec: next.psi, lane: next.lane, label: `${next.label}${NQ > 1 ? ` on q${next.lane}` : ""} (coming)` };
    return null;
  }

  function circuitChanged() {
    $("qiskitCode").textContent = qiskitCode();
    G.out = runCells(G.cells);
    updateMonitor();
    updateHUD();
  }

  // ---------- Panels ----------
  function buildTools() {
    const names = { X: "X", H: "H", Z: "Z", RX: "RX(θ)", RY: "RY(θ)", RZ: "RZ(θ)", CNOT: "CNOT", ERASE: "Eraser" };
    $("tools").innerHTML = [...lvl().tools, "ERASE"].map(t =>
      `<button data-tool="${t}" class="tool t-${t.toLowerCase()}${t === G.tool ? " on" : ""}">${names[t]}</button>`).join("");
    document.querySelectorAll("#tools button").forEach(b => b.addEventListener("click", () => {
      G.tool = b.dataset.tool;
      document.querySelectorAll("#tools button").forEach(x => x.classList.toggle("on", x === b));
      syncControls();
    }));
  }
  function buildBars() {
    const el = $("probBars");
    el.innerHTML = "";
    el.style.gridTemplateColumns = `repeat(${DIM}, 1fr)`;
    for (let i = 0; i < DIM; i++) {
      const bar = document.createElement("div");
      bar.className = "bar";
      bar.innerHTML = `<div class="track"><div class="fill"></div><span class="val"></span></div><span class="lbl">${i.toString(2).padStart(NQ, "0")}</span>`;
      el.appendChild(bar);
    }
  }
  function renderWave() {
    $("waveList").innerHTML = G.wave.map(w => `<span class="wv${w.type === "Q" ? " q" : ""}${w.done ? " done" : ""}">${w.type === "Q" ? "Q " : ""}${w.label}${NQ > 1 ? ` · q${w.lane}` : ""}</span>`).join("")
      + (G.boss || lvl().boss ? `<span class="wv boss">Boss ${G.boss ? G.boss.label : ""}</span>` : "");
  }
  function syncControls() {
    const rot = lvl().tools.some(t => ROT.includes(t));
    $("thetaRow").classList.toggle("hidden", !rot);
    $("thetaRow").classList.toggle("dim", !ROT.includes(G.tool));
    $("thetaVal").textContent = `θ = ${piFrac(G.theta)}`;
    $("theta").value = Math.round(G.theta / (PI / 12));
    $("toolHint").textContent = G.tool === "CNOT" ? "CNOT: click the control lane; the other lane is the target."
      : ROT.includes(G.tool) ? `${G.tool}(θ) uses the θ below. Click a planted rotation with a new θ to retune it.`
      : G.tool === "ERASE" ? "Click a planted gate to remove it." : `${G.tool}: click an open slot on a lane.`;
    document.querySelectorAll("#levelSel button").forEach((b, i) => b.classList.toggle("on", i === G.level));
    $("testLaneRow").classList.toggle("hidden", NQ < 2);
    document.querySelectorAll("#testLane button").forEach(b => b.classList.toggle("on", +b.dataset.l === G.testLane));
  }
  function updateMonitor() {
    const a = G.out || runCells(G.cells);
    $("outKet").textContent = fmtNamed(a);
    const p = probs(a);
    [...$("probBars").children].forEach((bar, i) => {
      bar.querySelector(".fill").style.height = `${p[i] * 100}%`;
      bar.querySelector(".val").textContent = `${Math.round(p[i] * 100)}%`;
    });
    const T = currentTarget();
    $("targetKet").textContent = T ? T.label : "No quantum target in this wave.";
    if (NQ === 2) {
      const c = concurrence(a);
      $("entInfo").classList.remove("hidden");
      $("entInfo").textContent = `Entanglement (concurrence) ${c.toFixed(2)}: ${c > 0.99 ? "maximally entangled" : c < 0.01 ? "product state, no entanglement" : "partly entangled"}. Each qubit's arrow shrinks toward the centre as entanglement grows.`;
    } else $("entInfo").classList.add("hidden");
    if (!G.showFid) $("fidInfo").textContent = "Match chance hidden. Turn it on to check your circuit, or find out in battle.";
    else if (!T) $("fidInfo").textContent = "No quantum target to compare against.";
    else {
      const f = T.kind === "boss" ? cabs2(inner(T.vec, a)) : overlapQubit(a, T.lane, T.vec);
      $("fidInfo").textContent = `Chance that one register passes the target's test: ${(f * 100).toFixed(1)}%.`
        + (T.kind === "boss" ? " The boss falls once your matches prove a rate above 50%; no un-entangled pair gets above 50%." : "");
    }
    drawBlochPanel(a, T);
    const r = [...G.regs].reverse().find(x => x.next > 0);
    $("regInfo").textContent = !r ? "No register on the lawn yet."
      : r.measured ? `Register #${r.id} measured: ${G.lastOutcome}.`
      : `Register #${r.id} at x = ${r.x.toFixed(2)}, ${r.next} of ${SLOTS.length} slots applied: ${fmtKet(r.amps)}`;
  }
  function updateHUD() {
    $("levelChip").textContent = `${lvl().key} · ${G.level + 1} of ${LEVELS.length}`;
    $("baseHP").textContent = `${G.baseHP} / 10`;
    $("gateHP").textContent = gateAlive() ? `${G.gateHP} / 10` : "destroyed";
    $("score").textContent = G.score;
    $("waveChip").textContent = `${G.waveIdx} / ${G.wave.length}${G.boss ? " + boss" : ""}`;
    $("pause").textContent = G.paused ? "Resume" : "Pause";
    $("basis").textContent = `Measurement basis: ${G.basis} (${G.basis === "Z" ? "0 / 1" : "+ / −"}), click to switch`;
    $("ready").classList.toggle("hidden", G.phase !== "plan");
    $("ready").textContent = `Ready: start level ${lvl().key}`;
    document.querySelectorAll("[data-z]").forEach(b => { b.disabled = G.phase !== "battle"; });
    $("entityCount").textContent = `Registers ${G.regs.length} · Zombies ${G.zombies.length}`;
  }

  // ---------- Bloch panel ----------
  function drawBlochPanel(a, T) {
    const c = bctx;
    c.clearRect(0, 0, 400, 200);
    const spheres = NQ === 1
      ? [{ cx: 200, title: "q0", vecs: [{ v: blochOf(a, 0), col: "#2f9e1f", name: "yours" },
        ...(T ? [{ v: blochOf(T.vec, 0), col: "#7a4fc0", name: "target" }] : [])] }]
      : lanes().map(q => ({ cx: 100 + q * 200, title: `q${q}`, vecs: [{ v: blochOf(a, q), col: "#2f9e1f", name: "yours" },
        ...(T && T.kind === "q" && T.lane === q ? [{ v: blochOf(T.vec, 0), col: "#7a4fc0", name: "target" }] : [])] }));
    for (const s of spheres) drawSphere(c, s.cx, 104, 70, s.vecs, s.title);
  }
  function drawSphere(c, cx, cy, R, vecs, title) {
    const P = v => [cx + R * (0.94 * v.y - 0.36 * v.x), cy - R * (0.94 * v.z - 0.22 * v.x)];
    const dark = matchMedia && matchMedia("(prefers-color-scheme: dark)").matches;
    const ink = dark ? "#f6e7c1" : "#3b2412", faint = dark ? "rgba(246,231,193,.55)" : "rgba(59,36,18,.4)";
    c.lineWidth = 1.5; c.strokeStyle = faint;
    c.beginPath(); c.arc(cx, cy, R, 0, PI * 2); c.stroke();
    c.beginPath();
    for (let k = 0; k <= 48; k++) { const t = (k / 48) * PI * 2; const [x, y] = P({ x: Math.cos(t), y: Math.sin(t), z: 0 }); k ? c.lineTo(x, y) : c.moveTo(x, y); }
    c.stroke();
    for (const [v, lab] of [[{ x: 0, y: 0, z: 1 }, "|0⟩"], [{ x: 0, y: 0, z: -1 }, "|1⟩"], [{ x: 1, y: 0, z: 0 }, "|+⟩"], [{ x: 0, y: 1, z: 0 }, "|+i⟩"]]) {
      const [x, y] = P(v); const [x2, y2] = P({ x: -v.x, y: -v.y, z: -v.z });
      c.beginPath(); c.moveTo(x2, y2); c.lineTo(x, y); c.stroke();
      const [lx, ly] = P({ x: v.x * 1.22, y: v.y * 1.22, z: v.z * 1.22 });
      c.font = `12px ${FONT}`; c.fillStyle = ink; c.textAlign = "center"; c.textBaseline = "middle";
      c.fillText(lab, lx, ly);
    }
    vecs.forEach((o, k) => {
      const [x, y] = P(o.v);
      c.strokeStyle = o.col; c.fillStyle = o.col; c.lineWidth = k ? 2.5 : 3.5;
      if (k) c.setLineDash([5, 4]);
      c.beginPath(); c.moveTo(cx, cy); c.lineTo(x, y); c.stroke(); c.setLineDash([]);
      c.beginPath(); c.arc(x, y, k ? 4 : 5, 0, PI * 2); c.fill();
    });
    c.font = `14px ${FONT}`; c.fillStyle = ink; c.textAlign = "left";
    c.fillText(title, cx - R - 20, cy - R + 4);
    if (vecs[0]) { c.font = `11px ${FONT}`; c.textAlign = "right"; c.fillText(`|r| = ${vecs[0].v.r.toFixed(2)}`, cx + R + 24, cy + R + 2); }
  }

  // ---------- Drawing ----------
  function rr(x, y, w, h, r) { ctx.beginPath(); ctx.roundRect(x, y, w, h, r); }
  function circle(x, y, r) { ctx.beginPath(); ctx.arc(x, y, r, 0, Math.PI * 2); }
  function text(s, x, y, size, color, align = "center") {
    ctx.font = `${size}px ${FONT}`; ctx.fillStyle = color; ctx.textAlign = align; ctx.textBaseline = "middle";
    ctx.fillText(s, x, y);
  }

  function drawLawn() {
    ctx.fillStyle = "#3d8a2a"; ctx.fillRect(0, 0, W, TOP);
    ctx.fillStyle = "#4c9c33";
    for (let x = 0; x < W; x += 34) { circle(x + 17, TOP - 2, 18); ctx.fill(); }
    for (let q = 0; q < NQ; q++) {
      for (let u = 0; u < LANE_END; u++) {
        ctx.fillStyle = (q + u) % 2 ? "#6cba43" : "#7cc952";
        ctx.fillRect(xPos(u), TOP + q * LANE_H, xPos(u + 1) - xPos(u) + 1, LANE_H);
      }
    }
    ctx.fillStyle = "rgba(255,255,240,.12)";
    ctx.fillRect(xPos(0.4), TOP, xPos(MEASURE_X) - xPos(0.4), DATA_H);
    const B = TOP + DATA_H;
    ctx.fillStyle = "#b98a4e"; ctx.fillRect(0, B, W, H - B);
    for (const u of [0, 5, 10, 15]) text(`x = ${u}`, xPos(u), H - 12, 12, "#5b3a1a");
    text("lab zone: quantum fights", (xPos(0.4) + xPos(MEASURE_X)) / 2, B + 17, 12, "#5b3a1a");
    text("beyond the gate: classical fights", (xPos(MEASURE_X) + xPos(LANE_END)) / 2, B + 17, 12, "#5b3a1a");
  }

  function drawHouse() {
    const hb = TOP + DATA_H;
    ctx.fillStyle = "#efdfb8"; ctx.fillRect(0, 0, 95, hb);
    ctx.strokeStyle = "#cdb98b"; ctx.lineWidth = 2;
    for (let y = 10; y < hb; y += 18) { ctx.beginPath(); ctx.moveTo(0, y); ctx.lineTo(95, y); ctx.stroke(); }
    ctx.fillStyle = "#8a5a2b"; ctx.fillRect(88, 0, 9, hb);
    const dy = TOP + DATA_H / 2;
    rr(24, dy - 38, 44, 76, 6); ctx.fillStyle = "#7a4a22"; ctx.fill();
    circle(60, dy, 4); ctx.fillStyle = "#ffd23f"; ctx.fill();
    text(`Base ${G.baseHP}`, 46, hb + 18, 15, "#fff4cf");
    if (G.houseFlash > 0) { ctx.fillStyle = `rgba(220,50,40,${G.houseFlash * 0.5})`; ctx.fillRect(0, 0, 97, hb); }
  }

  function drawSprout(x, y, t) {
    ctx.save(); ctx.translate(x, y + 6 + Math.sin(t * 2 + y) * 1.2); ctx.scale(0.58, 0.58);
    ctx.strokeStyle = "#2f7d1f"; ctx.lineWidth = 6;
    ctx.beginPath(); ctx.moveTo(0, 38); ctx.quadraticCurveTo(-6, 18, 0, 0); ctx.stroke();
    ctx.fillStyle = "#73d044"; ctx.strokeStyle = "#2f7d1f"; ctx.lineWidth = 3;
    rr(12, -26, 22, 20, 7); ctx.fill(); ctx.stroke();
    circle(0, -16, 20); ctx.fill(); ctx.stroke();
    circle(4, -24, 6); ctx.fillStyle = "#fff"; ctx.fill();
    circle(6, -24, 3); ctx.fillStyle = "#111"; ctx.fill();
    ctx.restore();
  }

  const TILE = { X: ["#ffb347", "#9a5b10", "#5b3a1a"], H: ["#dff6c8", "#2f7d1f", "#2f7d1f"], Z: ["#e9d7ff", "#5b3a8c", "#3b1f6e"],
    RX: ["#ffd6d0", "#a33a2a", "#7a2418"], RY: ["#d2f1ff", "#1f6f9a", "#124a6a"], RZ: ["#f3e3ff", "#6b3fa0", "#46206e"] };
  function gateTile(cx, y, g) {
    const [fill, stroke, ink] = TILE[g.type];
    const rot = ROT.includes(g.type), w = rot ? 44 : 30, h = rot ? 40 : 30;
    rr(cx - w / 2, y - h / 2, w, h, 6); ctx.fillStyle = fill; ctx.fill();
    ctx.lineWidth = 2.5; ctx.strokeStyle = stroke; ctx.stroke();
    if (rot) { text(g.type, cx, y - 8, 14, ink); text(piFrac(g.theta), cx, y + 9, 12, ink); }
    else text(g.type, cx, y + 1, 16, ink);
  }

  function drawCircuit(t) {
    const x1 = xPos(0.3), x2 = xPos(MEASURE_X);
    for (const q of lanes()) {
      ctx.strokeStyle = "rgba(60,40,15,.45)"; ctx.lineWidth = 2;
      ctx.beginPath(); ctx.moveTo(x1, laneY(q)); ctx.lineTo(x2, laneY(q)); ctx.stroke();
      drawSprout(xPos(0.17), laneY(q) - 4, t);
      text(`q${q}`, xPos(0.17) - 2, laneY(q) + 38, 12, "#2a4a12");
    }
    SLOTS.forEach(s => {
      const cx = xPos(s.u), cell = G.cells[s.i];
      if (cell.cx) {
        const yc = laneY(cell.cx.c), yt = laneY(cell.cx.t);
        ctx.strokeStyle = "#3b2412"; ctx.lineWidth = 3;
        ctx.beginPath(); ctx.moveTo(cx, yc); ctx.lineTo(cx, yt); ctx.stroke();
        circle(cx, yc, 8); ctx.fillStyle = "#3b2412"; ctx.fill();
        circle(cx, yt, 15); ctx.fillStyle = "#ffe9a6"; ctx.fill(); ctx.strokeStyle = "#3b2412"; ctx.lineWidth = 2.5; ctx.stroke();
        ctx.beginPath(); ctx.moveTo(cx - 15, yt); ctx.lineTo(cx + 15, yt); ctx.moveTo(cx, yt - 15); ctx.lineTo(cx, yt + 15); ctx.stroke();
        return;
      }
      for (const q of lanes()) {
        const g = cell.g[q];
        if (g) { gateTile(cx, laneY(q), g); continue; }
        ctx.save(); ctx.setLineDash([4, 4]); ctx.lineWidth = 1.5; ctx.strokeStyle = "rgba(255,255,255,.75)";
        rr(cx - 13, laneY(q) - 13, 26, 26, 5); ctx.stroke(); ctx.restore();
      }
    });
    const a = xPos(SLOTS[0].u) - 24, b = xPos(SLOTS[SLOTS.length - 1].u) + 24, y = 30;
    ctx.strokeStyle = "#c8402c"; ctx.lineWidth = 2;
    ctx.beginPath(); ctx.moveTo(a, y + 6); ctx.quadraticCurveTo(a, y, a + 8, y); ctx.lineTo(b - 8, y); ctx.quadraticCurveTo(b, y, b, y + 6); ctx.stroke();
    const s = "plant slots";
    ctx.font = `13px ${FONT}`;
    const tw = ctx.measureText(s).width + 14;
    rr((a + b) / 2 - tw / 2, y - 20, tw, 18, 5); ctx.fillStyle = "rgba(255,248,220,.94)"; ctx.fill();
    text(s, (a + b) / 2, y - 11, 13, "#c8402c");
  }

  function drawGate(t) {
    const x = xPos(MEASURE_X);
    if (!gateAlive()) {
      ctx.fillStyle = "#7a5a35";
      for (const q of lanes()) { ctx.save(); ctx.translate(x, laneY(q) + 30); ctx.rotate(0.5 - q * 0.4); ctx.fillRect(-12, -3, 24, 6); ctx.restore(); }
      text("gate destroyed", x, TOP + DATA_H + 16, 12, "#8a2a1a");
      return;
    }
    ctx.save();
    if (G.gateFlash > 0) { ctx.shadowColor = "#fff3a0"; ctx.shadowBlur = 50 * G.gateFlash; }
    rr(x - 12, TOP - 2, 24, DATA_H + 4, 6); ctx.fillStyle = "#9a6a35"; ctx.fill();
    ctx.restore();
    ctx.lineWidth = 3; ctx.strokeStyle = "#5b3a1a"; rr(x - 12, TOP - 2, 24, DATA_H + 4, 6); ctx.stroke();
    for (const q of lanes()) {
      const y = laneY(q);
      circle(x, y, 18); ctx.fillStyle = "#fff8dc"; ctx.fill(); ctx.strokeStyle = "#5b3a1a"; ctx.lineWidth = 2.5; ctx.stroke();
      const ang = -Math.PI / 2 + (G.gateFlash > 0 ? Math.sin(t * 40 + q) * 0.9 : 0.5 * Math.sin(t + q));
      ctx.beginPath(); ctx.moveTo(x, y + 5); ctx.lineTo(x + Math.cos(ang) * 12, y + 5 + Math.sin(ang) * 12);
      ctx.strokeStyle = "#c8402c"; ctx.lineWidth = 2.5; ctx.stroke();
    }
    rr(x - 52, TOP + DATA_H + 6, 104, 22, 6); ctx.fillStyle = "#8a5a2b"; ctx.fill();
    text(`Measure ${G.basis} · ${G.gateHP} HP`, x, TOP + DATA_H + 17, 12, "#fff4cf");
  }

  function drawRegs() {
    for (const r of G.regs) {
      const X = xPos(r.x);
      if (!r.measured && NQ > 1 && r.lanes.every(Boolean) && blochOf(r.amps, 0).r < 0.98) {
        ctx.save(); ctx.setLineDash([5, 4]); ctx.strokeStyle = "rgba(120,70,200,.85)"; ctx.lineWidth = 2.5;
        ctx.beginPath(); ctx.moveTo(X, laneY(0) - 8); ctx.lineTo(X, laneY(NQ - 1) - 8); ctx.stroke(); ctx.restore();
      }
      for (const q of lanes()) {
        if (!r.lanes[q]) continue;
        const y = laneY(q) - 8;
        ctx.save();
        if (r.measured) {
          circle(X, y, r.flash > 0 ? 10 : 8);
          ctx.fillStyle = r.basis === "Z" ? "#8be04e" : "#7fd8f0"; ctx.fill();
          ctx.lineWidth = 2; ctx.strokeStyle = "#2f5d1f"; ctx.stroke(); ctx.restore();
          text(sym(r.sym[q]), X, y - 17, 13, "#1d3d10");
        } else {
          ctx.shadowColor = "#b388ff"; ctx.shadowBlur = 12;
          circle(X, y, 8); ctx.fillStyle = "#d8c2ff"; ctx.fill();
          ctx.lineWidth = 2; ctx.strokeStyle = "#5b3a8c"; ctx.stroke(); ctx.restore();
          if (r.next >= SLOTS.length) text(shortLabel(r.amps, q), X, y - 17, 12, "#3b1f6e");
        }
      }
    }
  }

  function drawZombieBody(x, y, t, o) {
    const sw = o.still ? 0 : Math.sin(t * 2.6 + o.phase);
    ctx.save();
    ctx.translate(x, y); ctx.scale(o.scale, o.scale);
    const coat = o.flash ? "#ffffff" : o.coat, skin = o.flash ? "#ffffff" : o.skin;
    ctx.lineCap = "round";
    ctx.strokeStyle = "#3b3f5c"; ctx.lineWidth = 7;
    ctx.beginPath(); ctx.moveTo(-6, 0); ctx.lineTo(-9 + sw * 5, 28); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(6, 0); ctx.lineTo(8 - sw * 5, 28); ctx.stroke();
    rr(-15, -44, 30, 48, 8); ctx.fillStyle = coat; ctx.fill();
    ctx.strokeStyle = skin; ctx.lineWidth = 7;
    ctx.beginPath(); ctx.moveTo(-10, -34); ctx.lineTo(-36, -34 + sw * 3); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(-6, -26); ctx.lineTo(-32, -25 - sw * 3); ctx.stroke();
    circle(-2, -60, 17); ctx.fillStyle = skin; ctx.fill();
    circle(-9, -64, 5.5); ctx.fillStyle = "#fff"; ctx.fill(); circle(4, -63, 4.5); ctx.fill();
    ctx.fillStyle = "#222"; circle(-10, -63, 2.3); ctx.fill(); circle(3, -62, 2.1); ctx.fill();
    ctx.strokeStyle = "#3a2a1a"; ctx.lineWidth = 2;
    ctx.beginPath(); ctx.moveTo(-12, -51); ctx.lineTo(2, -52); ctx.stroke();
    if (o.crown) {
      ctx.fillStyle = "#ffd23f"; ctx.beginPath();
      ctx.moveTo(-16, -74); ctx.lineTo(-12, -92); ctx.lineTo(-4, -80); ctx.lineTo(2, -94); ctx.lineTo(8, -80); ctx.lineTo(14, -90); ctx.lineTo(12, -74);
      ctx.closePath(); ctx.fill(); ctx.strokeStyle = "#9a5b10"; ctx.stroke();
    }
    ctx.restore();
  }
  function card(s, x, y, bg) {
    ctx.font = `12px ${FONT}`;
    const w = Math.max(52, ctx.measureText(s).width + 14);
    const cx = Math.min(W - w / 2 - 4, x);
    rr(cx - w / 2, y - 9, w, 18, 5); ctx.fillStyle = bg; ctx.fill();
    text(s, cx, y, 12, "#fff4cf");
  }

  function drawZombies(t) {
    for (const z of G.zombies) {
      const x = xPos(z.x) - 8, y = laneY(z.lane) + 26, q = z.type === "Q";
      drawZombieBody(x, y, t, { scale: 0.8, phase: z.phase, flash: z.flash > 0, still: false,
        coat: q ? "#6b4fa0" : "#6e5a43", skin: q ? "#c9b8ef" : "#9fc27a" });
      card(q ? `${z.label} · ${"♥".repeat(Math.max(0, z.hp))}` : `${sym(z.type)}  ${z.hits}/2`, x, y - 84, q ? "#5b3a8c" : "#8a5a2b");
    }
  }

  function drawBoss(t) {
    const B = G.boss;
    if (!B) return;
    const x = xPos(B.x);
    ctx.save(); ctx.setLineDash([6, 5]); ctx.strokeStyle = "rgba(120,70,200,.9)"; ctx.lineWidth = 4;
    ctx.beginPath(); ctx.moveTo(x, laneY(0) - 30);
    for (let k = 1; k <= 20; k++) { const yy = laneY(0) - 30 + ((laneY(NQ - 1) - laneY(0)) * k) / 20; ctx.lineTo(x + Math.sin(t * 4 + k) * 5, yy); }
    ctx.stroke(); ctx.restore();
    for (const q of lanes()) {
      drawZombieBody(x, laneY(q) + 30, t, { scale: 1.0, phase: q * 1.7, flash: B.flash > 0,
        still: G.phase === "plan", coat: "#3a2a5e", skin: "#b9a2f0", crown: true });
    }
    // Card: the Bell state, then a proof bar. The bar is the 3σ lower bound on your match rate;
    // the boss falls when it passes the 50% mark (the best any un-entangled pair can do).
    const s = `Boss ${B.label}`, line = `matches ${B.matches} of ${B.tests} · prove more than 50%`;
    ctx.font = `14px ${FONT}`;
    const w1 = ctx.measureText(s).width;
    ctx.font = `12px ${FONT}`;
    const bw = Math.max(170, w1, ctx.measureText(line).width) + 24;
    const bx = Math.min(W - bw / 2 - 6, Math.max(bw / 2 + 6, x)), gx0 = bx - bw / 2 + 12, gw = bw - 24;
    rr(bx - bw / 2, 3, bw, 56, 7); ctx.fillStyle = "rgba(40,20,60,.92)"; ctx.fill();
    text(s, bx, 14, 14, "#efe6ff");
    rr(gx0, 26, gw, 8, 4); ctx.fillStyle = "#3a2a4a"; ctx.fill();
    rr(gx0, 26, gw * Math.max(0.02, wilsonLow(B.matches, B.tests)), 8, 4); ctx.fillStyle = "#b6f27a"; ctx.fill();
    ctx.fillStyle = "#ff6b5e"; ctx.fillRect(gx0 + gw * 0.5 - 1.5, 22, 3, 16);
    text(line, bx, 47, 12, "#efe6ff");
  }

  function drawSplats() {
    for (const s of G.splats) {
      ctx.save(); ctx.globalAlpha = Math.max(0, Math.min(1, s.life * 2.5));
      circle(s.x, s.y, s.r); ctx.strokeStyle = s.color; ctx.lineWidth = 3; ctx.stroke(); ctx.restore();
    }
  }

  // ---------- New-zombie callouts ----------
  // The first time each kind of zombie or boss shows up, a parchment tag with a red bracket (the
  // same style as the circuit tags) says what it is and what beats it. It follows the zombie for
  // CALLOUT_SECONDS of battle time, and stays up while you are still planning.
  const CALLOUT_SECONDS = 7;
  const CALLOUTS = {
    Z: { title: "Classical zombie (0 / 1)", lines: [
      "It fights beyond the measurement gate.",
      "Two measured bullets that read its symbol knock it out."] },
    X: { title: "Classical zombie (+ / −)", lines: [
      "Its symbol is + or −, so measure in the X basis.",
      "Two matching hits knock it out."] },
    Q: { title: "Quantum zombie", lines: [
      "It carries its own quantum state (on its card) in the lab zone.",
      "Only an unmeasured bullet in that same state can hit it.",
      "Classical bullets and other states pass straight through."] },
    boss: { title: "Entangled boss", lines: [
      "A Bell state spread over both lanes.",
      "Only an entangled pair of bullets can match it.",
      "Prove more than 50% matches: no unentangled pair can."] }
  };
  const zombieKind = z => (z.type === "Q" ? "Q" : z.type === "+" || z.type === "-" ? "X" : "Z");
  const bossKind = () => "boss";
  const zombieAnchor = z => ({ x: xPos(z.x) - 8, top: laneY(z.lane) - 67, bottom: laneY(z.lane) + 52 });
  const bossSpan = () => ({ x: xPos(G.boss.x) - 50, y1: laneY(0) - 72, y2: laneY(NQ - 1) + 56 });
  const seenKinds = new Set();
  let callouts = [], calloutClock = 0;
  function meetKind(kind, target, boss) {
    if (!kind || seenKinds.has(kind)) return;
    seenKinds.add(kind);
    callouts.push({ kind, target, boss, left: CALLOUT_SECONDS });
  }
  function tagSize(C) {
    ctx.font = `14px ${FONT}`;
    let w = ctx.measureText(`NEW · ${C.title}`).width;
    ctx.font = `12px ${FONT}`;
    for (const s of C.lines) w = Math.max(w, ctx.measureText(s).width);
    return { w: w + 20, h: 24 + 16 * C.lines.length };
  }
  function drawTag(C, x, y, m) {
    rr(x, y, m.w, m.h, 6); ctx.fillStyle = "rgba(255,248,220,.96)"; ctx.fill();
    ctx.strokeStyle = "rgba(200,64,44,.55)"; ctx.lineWidth = 1.5; ctx.stroke();
    text(`NEW · ${C.title}`, x + 10, y + 13, 14, "#c8402c", "left");
    C.lines.forEach((s, i) => text(s, x + 10, y + 30 + 16 * i, 12, "#5b3a1a", "left"));
  }
  // Tag above the zombie (below it when there is no room), with a bracket over its head.
  function zombieCallout(c) {
    if (!c.gone) c.a = zombieAnchor(c.target);
    const C = CALLOUTS[c.kind], a = c.a, m = tagSize(C);
    const above = a.top - 22 - m.h >= 4, d = above ? 6 : -6;
    const by = above ? a.top - 8 : a.bottom + 8, ty = above ? by - 14 - m.h : by + 14;
    const tx = Math.max(4, Math.min(W - m.w - 4, a.x - m.w / 2));
    ctx.strokeStyle = "#c8402c"; ctx.lineWidth = 2;
    ctx.beginPath(); ctx.moveTo(a.x - 26, by + d); ctx.quadraticCurveTo(a.x - 26, by, a.x - 18, by);
    ctx.lineTo(a.x + 18, by); ctx.quadraticCurveTo(a.x + 26, by, a.x + 26, by + d); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(a.x, by); ctx.lineTo(Math.max(tx + 12, Math.min(tx + m.w - 12, a.x)), above ? ty + m.h : ty); ctx.stroke();
    drawTag(C, tx, ty, m);
  }
  // Tag to the left of the boss, with a bracket along its whole height.
  function bossCallout(c) {
    if (!c.gone) c.s = bossSpan();
    const C = CALLOUTS[c.kind], s = c.s, m = tagSize(C), ym = (s.y1 + s.y2) / 2;
    const tx = Math.max(4, s.x - 16 - m.w), ty = Math.max(4, ym - m.h / 2);
    ctx.strokeStyle = "#c8402c"; ctx.lineWidth = 2;
    ctx.beginPath(); ctx.moveTo(s.x + 6, s.y1); ctx.quadraticCurveTo(s.x, s.y1, s.x, s.y1 + 8);
    ctx.lineTo(s.x, s.y2 - 8); ctx.quadraticCurveTo(s.x, s.y2, s.x + 6, s.y2); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(s.x, ym); ctx.lineTo(tx + m.w, ym); ctx.stroke();
    drawTag(C, tx, ty, m);
  }
  function drawCallouts(t) {
    const dt = Math.min(0.05, Math.max(0, t - (calloutClock || t)));
    calloutClock = t;
    for (const z of G.zombies || []) if (z.alive && !(z.wait > 0)) meetKind(zombieKind(z), z, false);
    if (G.boss) meetKind(bossKind(), G.boss, true);
    // A tag whose zombie has just fallen stays where it was for up to 3 more seconds.
    for (const c of callouts) {
      if (!c.gone && !(c.boss ? G.boss === c.target : c.target.alive && G.zombies.includes(c.target))) { c.gone = true; c.left = Math.min(c.left, 3); }
    }
    callouts = callouts.filter(c => c.left > 0);
    for (const c of callouts) {
      if ((G.phase === "battle" && !G.paused) || c.gone) c.left -= dt;
      ctx.save(); ctx.globalAlpha = Math.max(0, Math.min(1, c.left / 0.6));
      if (c.boss) bossCallout(c); else zombieCallout(c);
      ctx.restore();
    }
  }

  function draw(t) {
    ctx.clearRect(0, 0, W, H);
    drawLawn(); drawHouse(); drawCircuit(t); drawGate(t);
    drawZombies(t); drawBoss(t); drawRegs(); drawSplats(); drawCallouts(t);
    if (G.phase === "plan") {
      const a = xPos(6.0), b = xPos(12.2);
      rr(a, TOP + 10, b - a, 34, 8); ctx.fillStyle = "rgba(40,25,10,.8)"; ctx.fill();
      text(`Level ${lvl().key}: ${lvl().name}. Plant your lane, then press Ready.`, (a + b) / 2, TOP + 27, 14, "#fff4cf");
    }
    if (G.paused || G.phase === "over") {
      ctx.fillStyle = "rgba(20,40,10,.55)"; ctx.fillRect(0, 0, W, H);
      text(G.phase === "over" ? "Base lost" : "Paused", W / 2, H / 2 - 10, 40, "#fff4cf");
      text(G.phase === "over" ? "Restart to play again" : "Resume when ready", W / 2, H / 2 + 26, 16, "#fff4cf");
    }
  }

  // ---------- Events ----------
  $("theta").addEventListener("input", e => { G.theta = (+e.target.value * PI) / 12; syncControls(); });
  document.querySelectorAll("#levelSel button").forEach((b, i) => b.addEventListener("click", () => {
    if (G.phase === "battle" && !confirm("Leave the level in progress?")) return;
    loadLevel(i, false);
  }));
  document.querySelectorAll("#testLane button").forEach(b => b.addEventListener("click", () => { G.testLane = +b.dataset.l; syncControls(); }));
  document.querySelectorAll("[data-z]").forEach(b => b.addEventListener("click", () => {
    if (G.phase !== "battle") return;
    const type = b.dataset.z === "Q+" || b.dataset.z === "Qr" ? "Q" : b.dataset.z;
    const psi = type === "Q" ? stateFor(b.dataset.z === "Q+" ? "+" : "rand") : null;
    spawnZombie({ type, lane: G.testLane, psi, label: psi ? fmtNamed(psi) : sym(type) });
  }));
  $("clearCircuit").addEventListener("click", () => {
    G.cells = SLOTS.map(emptyCell);
    say("Lane cleared. Every register now reaches the gate as |" + "0".repeat(NQ) + "⟩."); circuitChanged();
  });
  $("basis").addEventListener("click", () => {
    if (G.phase === "over") return;
    G.basis = G.basis === "Z" ? "X" : "Z"; say(`Measurement basis switched to ${G.basis}.`); circuitChanged();
  });
  $("showFid").addEventListener("change", e => { G.showFid = e.target.checked; updateMonitor(); });
  $("ready").addEventListener("click", startWave);
  $("pause").addEventListener("click", () => { if (G.phase === "battle") { G.paused = !G.paused; updateHUD(); } });
  $("restart").addEventListener("click", () => loadLevel(0, false));

  canvas.addEventListener("click", ev => {
    const rc = canvas.getBoundingClientRect();
    const x = (ev.clientX - rc.left) * (W / rc.width), y = (ev.clientY - rc.top) * (H / rc.height);
    if (y < TOP || y >= TOP + DATA_H) return;
    const lane = Math.floor((y - TOP) / LANE_H);
    const s = SLOTS.find(ss => Math.abs(xPos(ss.u) - x) < 26);
    if (!s) return;
    const cell = G.cells[s.i];
    if (editCell(cell, lane)) {
      const g = cell.g[lane];
      say(cell.cx ? `CNOT planted in slot ${s.i + 1}: control q${cell.cx.c}, target q${cell.cx.t}.`
        : g ? `${gateText(g)} planted in slot ${s.i + 1}${NQ > 1 ? ` on q${lane}` : ""}.` : `Slot ${s.i + 1} cleared${NQ > 1 ? ` on q${lane}` : ""}.`);
      circuitChanged();
    }
  });

  // ---------- Tutorial storybook ----------
  // Pages shown before each level: a kid-friendly explanation, a "Try this" box and a collapsible
  // Physicist's corner. gt() draws a gate the way it looks on the lawn.
  const gt = (s, k) => `<span class="gt gt-${k}">${s}</span>`;
  const TUT = {
    intro: [{
      title: "Welcome to the Quantum Garden!",
      kid: `<p>Zombies are marching across your lawn toward the house. Your sprouts don't shoot peas. They shoot <b>qubits</b>: tiny pieces of quantum stuff.</p>
        <p>Before a qubit reaches the wooden <b>measurement gate</b> in the middle of the lawn, you can change it with the <b>gates</b> you plant in the slots. Change it the right way and it knocks the zombies out.</p>`,
      try: ["Read each page, then press <b>Let's play!</b>", "You can open these pages again at any time with <b>How to play</b>."],
      phys: `<p>Every register is an exact statevector simulated in JavaScript. Each planted gate is a unitary applied as the register crosses its slot, and the measurement gate is a projective measurement sampled with the Born rule. Every outcome comes from |amplitude|², not from a random-number shortcut.</p>`
    }],
    levels: [
      [{
        title: "A qubit is like a spinning coin",
        kid: `<p>A normal computer bit is like a light switch: it is either <b>0</b> or <b>1</b>.</p>
          <p>A qubit is more like a <b>spinning coin</b>. While it spins it isn't heads or tails yet; it is a bit of both at once. That is called <b>superposition</b>.</p>
          <p>When you <b>measure</b> it, it's like slapping the coin flat on the table: it has to choose 0 or 1, and it stays that way.</p>`,
        phys: `<p>A qubit state is |ψ⟩ = α|0⟩ + β|1⟩ with |α|² + |β|² = 1. Measuring in the Z basis gives 0 with probability |α|² and 1 with probability |β|² (the Born rule), and the state collapses onto the outcome.</p>`
      }, {
        title: "Your three gates",
        kid: `<p>${gt("X", "x")} <b>flips</b> the qubit: 0 becomes 1 and 1 becomes 0.</p>
          <p>${gt("H", "h")} <b>starts the coin spinning</b>. From 0 it makes a perfect half-and-half state called <b class="bk">|+⟩</b>.</p>
          <p>${gt("Z", "z")} is sneaky. It doesn't change the chances of 0 and 1 at all, but it flips the hidden <b>direction of the spin</b>, turning |+⟩ into <b class="bk">|−⟩</b>.</p>`,
        try: [`Plant ${gt("H", "h")}: the Qubit monitor shows |+⟩.`, `Plant ${gt("H", "h")} then ${gt("Z", "z")}: now it shows |−⟩, but the 0 and 1 bars stay at 50% each.`],
        phys: `<p>X = [[0, 1], [1, 0]], H = (1/√2)[[1, 1], [1, −1]], Z = diag(1, −1). |±⟩ = (|0⟩ ± |1⟩)/√2. Z|+⟩ = |−⟩ changes only the relative phase, which a Z measurement cannot see and an X measurement can.</p>`
      }, {
        title: "The measurement gate asks one question",
        kid: `<p>The wooden gate ${gt("Measure", "m")} asks every qubit a question.</p>
          <p>In <b>Z mode</b> it asks "0 or 1?". In <b>X mode</b> it asks "+ or −?".</p>
          <p>After the gate, qubits turn into ordinary <b>peas</b> with the answer printed on them: green peas for Z answers, blue peas for X answers. If zombies smash the gate (it has 10 HP), qubits stay quantum all the way down the lane.</p>`,
        try: ["Click the <b>Measurement basis</b> button to switch between Z and X."],
        phys: `<p>An X-basis measurement is H followed by a Z measurement: |+⟩ always gives +, |−⟩ always gives −. Measuring |+⟩ in Z gives a fair coin, while measuring it in X gives + every time. Choosing the basis is choosing which question you ask.</p>`
      }, {
        title: "Two kinds of zombies",
        kid: `<p><b>Classical zombies</b> (brown) wear a label: 0, 1, + or −. They live past the gate. Hit one <b>twice</b> with peas that carry the same label and it falls. A classical zombie that reaches the gate breaks itself on it, and the gate loses 1 HP.</p>
          <p><b>Quantum zombies</b> (purple) carry a quantum state on their card, like |+⟩. They walk straight through the gate, so only <b>unmeasured</b> qubits can fight them. If your qubit matches the zombie's state, the zombie loses a heart. If the states are different, they simply don't interact: your qubit passes straight through.</p>`,
        phys: `<p>A quantum zombie with state |z⟩ performs the projective test {|z⟩⟨z|, 1 − |z⟩⟨z|} on the qubit it meets. It matches with probability |⟨z|ψ⟩|²; two matches kill it. On a mismatch the qubit continues in the post-measurement state (1 − |z⟩⟨z|)|ψ⟩, normalised. An orthogonal qubit is left exactly unchanged: states that are fully different never interact.</p>`
      }, {
        title: "Your first battle",
        kid: `<p>Click a gate in the <b>Gate tray</b>, then click an empty slot on the lane. Press <b>Ready</b> when you're set. You can change gates and the measurement mode while the wave is moving.</p>
          <p>The <b>Incoming</b> bar shows which zombies are on their way.</p>`,
        try: ["Zombie <b>0</b> coming? No gates and Z mode: every pea is a 0.",
          `Zombie <b>1</b>? Plant ${gt("X", "x")}.`,
          `Zombie <b>+</b> or <b>−</b>? Plant ${gt("H", "h")} (add ${gt("Z", "z")} for −) and switch to X mode.`,
          `Quantum zombie <b>|+⟩</b>? Plant ${gt("H", "h")}: your unmeasured qubits match it perfectly.`,
          `Clever trick: ${gt("H", "h")} with Z mode fires random 0s and 1s, so one lane can hit both 0 and 1 zombies.`],
        phys: `<p>Your planted lane is a fixed one-qubit circuit and each register is an independent run of it. The strategy is to choose a state and a measurement basis whose outcome distribution matches the incoming wave.</p>`
      }],
      [{
        title: "The qubit globe",
        kid: `<p>Picture a qubit as an <b>arrow inside a globe</b>. The north pole is <b class="bk">|0⟩</b>, the south pole is <b class="bk">|1⟩</b>, and every point on the equator is a perfect half-and-half state, like |+⟩ or |−⟩.</p>
          <p>The globe in the <b>Qubit monitor</b> shows your arrow in <b>green</b> and the next quantum zombie's arrow in <b>purple</b>. Make the two arrows point the same way and you win!</p>`,
        phys: `<p>|ψ⟩ = cos(θ/2)|0⟩ + e^(iφ) sin(θ/2)|1⟩ corresponds to the Bloch vector (sin θ cos φ, sin θ sin φ, cos θ). A global phase does not move the arrow, so it never matters.</p>`
      }, {
        title: "Turning dials: rotation gates",
        kid: `<p>${gt("RY(θ)", "r")} tips the arrow from the north pole toward the south pole.</p>
          <p>${gt("RZ(θ)", "r")} spins the arrow around the globe, like a merry-go-round.</p>
          <p>${gt("RX(θ)", "r")} tips it too, but in a different direction.</p>
          <p>Choose θ with the <b>Angle</b> slider. It moves in steps of π/12, which is 15°.</p>`,
        phys: `<p>RY(θ)|0⟩ = cos(θ/2)|0⟩ + sin(θ/2)|1⟩. RZ(φ) = diag(e^(−iφ/2), e^(iφ/2)) adds the relative phase e^(iφ). RX(θ)|0⟩ = cos(θ/2)|0⟩ − i sin(θ/2)|1⟩.</p>`
      }, {
        title: "Reading a zombie's card",
        kid: `<p>A card like <span class="bk">√3/2|0⟩ + 1/2|1⟩</span> is a recipe. Square each number to get the chances: (√3/2)² = 3/4 for 0 and (1/2)² = 1/4 for 1.</p>
          <p>To build it, use ${gt("RY(θ)", "r")}: it puts cos(θ/2) in front of |0⟩. cos(30°) = √3/2, so θ/2 = 30° and θ = <b>π/3</b>.</p>
          <p>An <b>i</b> in front of |1⟩ means a quarter-turn spin: add ${gt("RZ(π/2)", "r")}. A <b>−</b> means a half turn: ${gt("RZ(π)", "r")}. A <b>−i</b> means ${gt("RZ(3π/2)", "r")}.</p>`,
        try: [`${gt("RY(π/3)", "r")} then ${gt("RZ(π/2)", "r")} makes √3/2|0⟩ + i/2|1⟩.`,
          "Tick <b>Show match chance</b>: 100% means a perfect copy of the zombie's state.",
          "A quantum zombie has 2 hearts. A mismatched qubit just flies past it, so the closer your state, the faster it falls."],
        phys: `<p>Each encounter passes with probability F = |⟨z|ψ⟩|² (the fidelity). Mismatched qubits pass through without effect, so the expected number of registers needed to remove both hearts is 2/F: a perfect copy (F = 1) needs exactly 2, and an orthogonal state (F = 0) never wins.</p>`
      }],
      [{
        title: "Magic twin coins",
        kid: `<p>Now you have <b>two lanes</b>, q0 and q1. Two qubits can become <b>entangled</b>. That's like two magic coins that always land the same way, even though each one on its own is completely random.</p>
          <p>You can't describe one twin without the other. On the lawn, a dashed purple line links entangled qubits. In the Qubit monitor, each qubit's arrow <b>shrinks toward the centre</b> of its globe.</p>`,
        phys: `<p>The Bell state |Φ+⟩ = (|00⟩ + |11⟩)/√2 cannot be written as a product |a⟩⊗|b⟩. Each qubit's reduced state is maximally mixed, so its Bloch vector has length 0. The concurrence C = 2|α₀₀α₁₁ − α₀₁α₁₀| equals 1.</p>`
      }, {
        title: "How to make twins",
        kid: `<p>${gt("CNOT", "c")} means: "if the <b>control</b> qubit is 1, flip the <b>target</b>." Click the lane you want as control; the other lane becomes the target.</p>
          <p>Recipe: ${gt("H", "h")} on q0, then ${gt("CNOT", "c")} with q0 as control. That makes the twins <b class="bk">(|00⟩ + |11⟩)/√2</b>.</p>`,
        try: [`Φ+ = (|00⟩ + |11⟩)/√2: ${gt("H", "h")} on q0, then ${gt("CNOT", "c")}.`,
          `Φ− = (|00⟩ − |11⟩)/√2: then add ${gt("Z", "z")} on q0.`,
          `Ψ+ = (|01⟩ + |10⟩)/√2: put ${gt("X", "x")} on q1 in the same slot as the H.`,
          `Ψ− = (|01⟩ − |10⟩)/√2: ${gt("X", "x")} on q1 first, and ${gt("Z", "z")} on q0 at the end.`],
        phys: `<p>CNOT(H⊗I)|00⟩ = |Φ+⟩. The four Bell states differ only by local X and Z gates, so one entangling gate plus single-qubit gates reaches all of them.</p>`
      }, {
        title: "The entangled boss",
        kid: `<p>The boss is a pair of crowned zombies tied together by a purple string. Its card shows a twin state such as (|00⟩ + |11⟩)/√2.</p>
          <p>It pushes through the gate and then tests <b>whole, unmeasured</b> registers in the lab zone, once they have passed all your gates. A matching register is used up; a different one passes straight through. The boss never heals.</p>
          <p>Here's the twist: a pair that is <b>not</b> entangled can match at most half the time. The boss only falls when your matches <b>prove</b> you beat one half. The green bar on its card shows how sure that proof is, and the red mark is the 50% line. <b>Only twins can beat this boss!</b></p>`,
        try: ["Build the twin state shown on the boss's card, in the first slots so it is ready early.",
          "Perfect twins prove it after 10 matches in a row.",
          "If the boss reaches the house it costs 3 HP and walks back, but your matches so far still count."],
        phys: `<p>For any product state, |⟨Φ|a⊗b⟩|² ≤ 1/2: the largest squared Schmidt coefficient of a maximally entangled state. So a Bell-state fidelity above 1/2 is an entanglement witness. The boss falls when the 3σ Wilson lower bound on the observed match rate exceeds 1/2; for a perfect Bell state (all matches) that takes n/(n + 9) > 1/2, i.e. 10 matches.</p>`
      }]
    ]
  };
  const tut = { pages: [], i: 0, seen: new Set(), resume: false };
  function openTutorial(withIntro) {
    const tag = `Level ${lvl().key}`;
    tut.pages = [...(withIntro ? TUT.intro.map(p => ({ ...p, tag: "Welcome" })) : []),
      ...(TUT.levels[G.level] || []).map(p => ({ ...p, tag }))];
    if (!tut.pages.length) return;
    tut.i = 0;
    tut.resume = G.phase === "battle" && !G.paused;
    if (tut.resume) G.paused = true;
    renderTutorial();
    $("tut").classList.remove("hidden");
    $("tutNext").focus();
  }
  function closeTutorial() {
    $("tut").classList.add("hidden");
    if (tut.resume) { G.paused = false; tut.resume = false; }
    updateHUD();
  }
  function renderTutorial() {
    const p = tut.pages[tut.i], n = tut.pages.length;
    $("tutStep").textContent = `${p.tag} · page ${tut.i + 1} of ${n}`;
    $("tutTitle").textContent = p.title;
    $("tutBody").innerHTML = p.kid + (p.try ? `<div class="tut-try"><b>Try this</b><ul>${p.try.map(t => `<li>${t}</li>`).join("")}</ul></div>` : "");
    $("tutPhysWrap").classList.toggle("hidden", !p.phys);
    $("tutPhys").innerHTML = p.phys || "";
    $("tutPhysWrap").open = false;
    $("tutBack").disabled = tut.i === 0;
    $("tutNext").textContent = tut.i === n - 1 ? "Let's play!" : "Next";
    $("tutDots").innerHTML = tut.pages.map((_, k) => `<i class="${k === tut.i ? "on" : ""}"></i>`).join("");
  }
  // Each level's pages open once per session, the first time it loads (with the welcome page before level 1).
  function autoTutorial() {
    if (tut.seen.has(G.level)) return;
    const first = tut.seen.size === 0;
    tut.seen.add(G.level);
    openTutorial(first);
  }
  $("tutNext").addEventListener("click", () => { if (tut.i < tut.pages.length - 1) { tut.i++; renderTutorial(); } else closeTutorial(); });
  $("tutBack").addEventListener("click", () => { if (tut.i > 0) { tut.i--; renderTutorial(); } });
  $("tutSkip").addEventListener("click", closeTutorial);
  $("tutOpen").addEventListener("click", () => openTutorial(true));
  document.addEventListener("keydown", e => {
    if ($("tut").classList.contains("hidden")) return;
    if (e.key === "Escape") closeTutorial();
    else if (e.key === "ArrowRight") $("tutNext").click();
    else if (e.key === "ArrowLeft") $("tutBack").click();
  });

  // ---------- Loop ----------
  let last = 0;
  function frame(ms) {
    const t = ms / 1000, dt = Math.min(t - (last || t), 0.05);
    last = t;
    update(dt);
    draw(t);
    G.monT += dt;
    if (G.monT > 0.15) { G.monT = 0; updateMonitor(); updateHUD(); }
    requestAnimationFrame(frame);
  }

  window.QPVZ1 = { G, update, loadLevel, startWave, editCell, circuitChanged, runCells, fmtKet, fmtNamed, BELL, qubitState };
  window.addEventListener("resize", fitCanvas);
  loadLevel(0, false);
  requestAnimationFrame(frame);
})();
