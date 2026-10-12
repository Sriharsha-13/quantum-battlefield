// Quantum Plants vs. Zombies: Adventure 2, "Grover's search" (4 levels)
// a. 2 data qubits + ancilla: build the oracle, diffuser locked. Kill bar > 50%.
// b. Build oracle and diffuser yourself. Kill bar > 70%.
// c. Superposed boss: 2 marked patterns out of 4. Standard Grover is stuck at exactly 50% for any
//    round count; phase-matched Grover (Long 2001) with phase i in oracle and diffuser reaches 100%.
//    Phase i comes from a controlled-√X onto the |−⟩ ancilla (√X|−⟩ = i|−⟩) and a controlled-S.
// d. 3 data qubits + ancilla, 1 of 8. Kill bar > 90%, so one round (78.1%) is not enough: repeat.
// Lane j is the phase-kickback ancilla, prepared in |−⟩ by locked X·H before the Grover loop and
// returned to |1⟩ by a locked H after it. Bit strings read q(n-1) ... q0, as in Qiskit.
(() => {
  "use strict";

  // ==QUANTUM CORE==
  let ND = 2, ANC = 2, DIM = 8, DDIM = 4;
  const setDims = nd => { ND = nd; ANC = nd; DIM = 1 << (nd + 1); DDIM = 1 << nd; };
  const SQ = Math.SQRT1_2;
  const C = (re, im = 0) => [re, im];
  const cadd = (a, b) => [a[0] + b[0], a[1] + b[1]];
  const cmul = (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]];
  const cconj = a => [a[0], -a[1]];
  const cabs2 = a => a[0] * a[0] + a[1] * a[1];
  const basisState = i => { const a = Array.from({ length: DIM }, () => C(0)); a[i] = C(1); return a; };
  const bitOf = (i, q) => (i >> q) & 1;
  const label = d => d.toString(2).padStart(ND, "0");
  const gX = (a, q) => a.map((_, i) => a[i ^ (1 << q)].slice());
  const gH = (a, q) => {
    const b = a.map(v => v.slice());
    for (let i = 0; i < DIM; i++) {
      if (bitOf(i, q)) continue;
      const j = i | (1 << q);
      b[i] = [(a[i][0] + a[j][0]) * SQ, (a[i][1] + a[j][1]) * SQ];
      b[j] = [(a[i][0] - a[j][0]) * SQ, (a[i][1] - a[j][1]) * SQ];
    }
    return b;
  };
  const gHdata = a => { let b = a; for (let q = 0; q < ND; q++) b = gH(b, q); return b; };
  const allSet = (i, qs) => qs.every(q => bitOf(i, q));
  // ⊕ onto j (CNOT / Toffoli / C3X): with j in |−⟩ this kicks back a phase −1.
  const gMCX = (a, qs) => a.map((v, i) => (allSet(i, qs) ? a[i ^ (1 << ANC)] : v).slice());
  // Controlled-Z family: phase −1 on all-ones of the chosen qubits.
  const gMCZ = (a, qs) => a.map((v, i) => (allSet(i, qs) ? [-v[0], -v[1]] : v.slice()));
  // Controlled-S family: phase i on all-ones of the chosen qubits.
  const gMCS = (a, qs) => a.map((v, i) => (allSet(i, qs) ? [-v[1], v[0]] : v.slice()));
  // Controlled-√X onto j: √X = ½[[1+i, 1−i], [1−i, 1+i]] and √X|−⟩ = i|−⟩, so this kicks back phase i.
  function gMCSX(a, qs) {
    const b = a.map(v => v.slice()), p = C(0.5, 0.5), m = C(0.5, -0.5);
    for (let i = 0; i < DIM; i++) {
      if (!allSet(i, qs) || bitOf(i, ANC)) continue;
      const j = i | (1 << ANC);
      b[i] = cadd(cmul(p, a[i]), cmul(m, a[j]));
      b[j] = cadd(cmul(m, a[i]), cmul(p, a[j]));
    }
    return b;
  }
  const CTRL = { MCX: gMCX, MCZ: gMCZ, MCS: gMCS, MCSX: gMCSX };

  function applyCell(a, cell) {
    if (!cell) return a;
    if (cell.ctrl) return CTRL[cell.ctrl.type](a, cell.ctrl.qs);
    let b = a;
    cell.x.forEach((on, q) => { if (on) b = gX(b, q); });
    return b;
  }
  const prep = a => gH(gX(gHdata(a), ANC), ANC);
  function grRound(a, oracle, diff) {
    let b = a;
    for (const c of oracle) b = applyCell(b, c);
    b = gHdata(b);
    for (const c of diff) b = applyCell(b, c);
    return gHdata(b);
  }
  function fullCircuit(oracle, diff, rounds) {
    let a = prep(basisState(0));
    for (let r = 0; r < rounds; r++) a = grRound(a, oracle, diff);
    return gH(a, ANC);
  }
  const dataProbs = a => Array.from({ length: DDIM }, (_, d) => cabs2(a[d]) + cabs2(a[d | (1 << ANC)]));
  function sampleIndex(p) {
    let r = Math.random() * p.reduce((s, v) => s + v, 0);
    for (let i = 0; i < p.length; i++) { r -= p[i]; if (r < 0) return i; }
    return p.length - 1;
  }
  function measureData(a, basis) {
    const b = basis === "X" ? gHdata(a) : a;
    const d = sampleIndex(dataProbs(b));
    const AB = 1 << ANC;
    const n = Math.sqrt(cabs2(b[d]) + cabs2(b[d | AB])) || 1;
    const c = Array.from({ length: DIM }, () => C(0));
    c[d] = [b[d][0] / n, b[d][1] / n]; c[d | AB] = [b[d | AB][0] / n, b[d | AB][1] / n];
    return { d, amps: basis === "X" ? gHdata(c) : c };
  }
  // Reduced Bloch vector of qubit q; its length drops below 1 when q is entangled.
  function bloch(a, q) {
    let rz = 0, r01 = C(0);
    for (let i = 0; i < DIM; i++) {
      rz += bitOf(i, q) ? -cabs2(a[i]) : cabs2(a[i]);
      if (!bitOf(i, q)) r01 = cadd(r01, cmul(a[i], cconj(a[i | (1 << q)])));
    }
    const rx = 2 * r01[0], ry = -2 * r01[1];
    return { rz, rx, ry, r: Math.hypot(rz, rx, ry) };
  }
  function measureQubitX(a, q) {
    const b = gH(a, q);
    let p0 = 0;
    for (let i = 0; i < DIM; i++) if (!bitOf(i, q)) p0 += cabs2(b[i]);
    const plus = Math.random() < p0;
    const c = b.map((v, i) => (bitOf(i, q) === (plus ? 0 : 1) ? v.slice() : C(0)));
    const n = Math.sqrt(c.reduce((s, v) => s + cabs2(v), 0)) || 1;
    return { plus, amps: gH(c.map(v => [v[0] / n, v[1] / n]), q) };
  }
  // ==END CORE==

  // threshold: the boss dies only if more than this fraction of shots hit its pattern(s).
  // diff: "locked" plants the textbook diffuser and locks it; "open" leaves it to the player.
  const LEVELS = [
    { nd: 2, key: "a", name: "Mark the boss", escorts: [], diff: "locked", pair: false, threshold: 0.5,
      tools: ["X", "MCX", "ERASE"], textbook: false,
      intro: "Two data qubits, four attack patterns. The diffuser is planted and locked; build the oracle that marks the boss. Kill bar: more than 500 of 1000 shots." },
    { nd: 2, key: "b", name: "Build the diffuser", escorts: ["0", "1"], diff: "open", pair: false, threshold: 0.7,
      tools: ["X", "MCX", "MCZ", "ERASE"], textbook: false,
      intro: "Now plant the diffuser too. Kill bar: more than 700 of 1000 shots. Try different repeat counts: with four patterns, more rounds do not mean more bullets." },
    { nd: 2, key: "c", name: "Superposed boss", escorts: ["Q", "+"], diff: "open", pair: true, threshold: 0.5,
      tools: ["X", "MCX", "MCZ", "MCSX", "MCS", "ERASE"], textbook: true,
      intro: "This boss hides in two patterns at once. Standard Grover stays stuck at exactly 50% for every repeat count, which is not enough. Find the phase that fixes it." },
    { nd: 3, key: "d", name: "Three lanes", escorts: ["0", "Q", "1"], diff: "open", pair: false, threshold: 0.9,
      tools: ["X", "MCX", "MCZ", "ERASE"], textbook: true,
      intro: "Three data qubits, eight patterns, a 4-qubit circuit with the ancilla. Kill bar: more than 900 of 1000 shots, so one Grover round (78%) is not enough." }
  ];
  const TOOL_HINT = {
    X: "X: click open slots on data lanes. A column holds X gates or one controlled gate.",
    MCX: "⊕ onto j. With j in |−⟩ it flips the sign of patterns whose controls are all 1. 1 control: CNOT · 2: Toffoli · 3: C3X.",
    MCZ: "Controlled-Z on the chosen qubits: phase −1 when all are 1. 1: Z · 2: CZ · 3: CCZ.",
    MCSX: "√X onto j. √X|−⟩ = i|−⟩, so this kicks back phase i instead of −1. 1 control: C-√X · 2: CC-√X.",
    MCS: "Controlled-S on the chosen qubits: phase i when all are 1. 1: S · 2: CS · 3: CCS.",
    ERASE: "Eraser: click a planted gate to remove it."
  };

  const $ = id => document.getElementById(id);
  const canvas = $("lawn");
  const ctx = canvas.getContext("2d");

  const W = 1200, H = 520, TOP = 44, DATA_H = 330, GAP = 14, POOL_H = 66;
  const POOL_TOP = TOP + DATA_H + GAP, BOTTOM = POOL_TOP + POOL_H;
  const X0 = 100, XM = 584, XE = 1172;
  const FONT = '"Lilita One", "Arial Rounded MT Bold", "Trebuchet MS", system-ui, sans-serif';
  const LANE_END = 15, MEASURE_X = 5, QSPEED = 1, EMIT = 2, ZSPEED = 0.45, BOSS_SPEED = 0.2, RAD = 0.3;
  let LANE_H = DATA_H / ND;
  let COLS = [];

  const xPos = u => (u <= MEASURE_X ? X0 + (u * (XM - X0)) / MEASURE_X
    : XM + ((u - MEASURE_X) * (XE - XM)) / (LANE_END - MEASURE_X));
  const laneY = q => (q === ANC ? POOL_TOP + POOL_H / 2 : TOP + q * LANE_H + LANE_H / 2);
  const emptyCell = () => ({ x: new Array(ND).fill(false), ctrl: null });
  const dataQs = () => Array.from({ length: ND }, (_, q) => q);

  function buildCols() {
    const hData = Object.fromEntries(dataQs().map(q => [q, "H"]));
    COLS = [
      { kind: "L", tiles: { ...hData, [ANC]: "X" }, apply: a => gX(gHdata(a), ANC) },
      { kind: "L", tiles: { [ANC]: "H" }, apply: a => gH(a, ANC) },
      { kind: "P", reg: "O", i: 0 }, { kind: "P", reg: "O", i: 1 }, { kind: "P", reg: "O", i: 2 },
      { kind: "L", tiles: hData, apply: gHdata },
      { kind: "P", reg: "D", i: 0 }, { kind: "P", reg: "D", i: 1 }, { kind: "P", reg: "D", i: 2 },
      { kind: "L", tiles: hData, apply: gHdata, blockEnd: true },
      { kind: "L", tiles: { [ANC]: "H" }, apply: a => gH(a, ANC) }
    ].map((c, i) => ({ ...c, u: 0.65 + 0.4 * i }));
  }

  const G = {
    level: 0, oracle: [], diff: [], rounds: 1, tool: "X", ctrlQs: [], lane: 0, basis: "Z", showPred: false,
    phase: "plan", baseHP: 10, gateHP: 10, score: 0,
    regs: [], zombies: [], boss: null, splats: [],
    emitT: 0, paused: false, nextId: 1, t: 0, loopFlash: 0, gateFlash: 0, houseFlash: 0, monT: 0,
    lastOutcome: ""
  };
  const lvl = () => LEVELS[G.level];
  const cellOf = c => (c.reg === "O" ? G.oracle[c.i] : G.diff[c.i]);
  const gateAlive = () => G.gateHP > 0;
  const ket = d => `|${label(d)}⟩`;
  const bossKet = B => (B.patterns.length === 1 ? ket(B.patterns[0]) : `(${B.patterns.map(ket).join(" + ")})/√2`);
  // Lower end of the 3σ (Wilson) confidence interval for a hit rate, from k hits in n shots.
  function wilsonLow(k, n, z = 3) {
    if (!n) return 0;
    const p = k / n, z2 = z * z;
    return (p + z2 / (2 * n) - z * Math.sqrt(p * (1 - p) / n + z2 / (4 * n * n))) / (1 + z2 / n);
  }
  const pct = t => `${Math.round(t * 100)}%`;
  const diffLocked = () => lvl().diff === "locked";
  const sym = s => (s === "-" ? "−" : s);
  const textbookDiff = () => [
    { x: new Array(ND).fill(true), ctrl: null },
    { x: new Array(ND).fill(false), ctrl: { type: "MCZ", qs: dataQs() } },
    { x: new Array(ND).fill(true), ctrl: null }];

  function fitCanvas() {
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    canvas.width = W * dpr; canvas.height = H * dpr;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  }
  function say(t) { $("msg").textContent = t; }

  // ---------- Levels ----------
  function loadLevel(i, keepStats) {
    G.level = i;
    const L = lvl();
    setDims(L.nd);
    LANE_H = DATA_H / ND;
    buildCols();
    G.oracle = [emptyCell(), emptyCell(), emptyCell()];
    G.diff = diffLocked() ? textbookDiff() : [emptyCell(), emptyCell(), emptyCell()];
    if (!L.tools.includes(G.tool)) G.tool = "X";
    G.ctrlQs = dataQs();
    G.rounds = 1;
    G.lane = 0;
    if (!keepStats) { G.baseHP = 10; G.score = 0; }
    G.gateHP = 10;
    G.paused = false;
    G.splats = [];
    $("duelInfo").textContent = "No duel yet.";
    ["ampBars", "predBars", "duelBars"].forEach(id => buildBars($(id), id === "ampBars"));
    renderBars($("duelBars"), new Array(DDIM).fill(0), new Array(DDIM).fill(""), []);
    syncControls();
    previewLevel();
    autoTutorial();
  }

  function previewLevel() {
    const L = lvl();
    G.phase = "plan";
    G.regs = [];
    G.zombies = [];
    let patterns;
    if (L.pair) {
      // Two patterns that agree on one data bit: the boss is spread evenly over the other bit.
      const q = Math.floor(Math.random() * ND), v = Math.floor(Math.random() * 2);
      patterns = Array.from({ length: DDIM }, (_, d) => d).filter(d => bitOf(d, q) === v);
    } else patterns = [1 + Math.floor(Math.random() * (DDIM - 1))];
    G.boss = { patterns, x: 13.6, prevX: 13.6, hits: 0, tests: 0, flash: 0, hitGate: false };
    L.escorts.forEach((type, k) => {
      G.zombies.push({ id: G.nextId++, type, lane: Math.floor(Math.random() * ND), x: LANE_END, prevX: LANE_END,
        hits: 0, alive: true, flash: 0, phase: Math.random() * 6, wait: 3 + k * 4 });
    });
    say(`Level ${L.key}: ${L.name}. ${L.intro} Boss ${bossKet(G.boss)}. Press Ready when your circuit is planted.`);
    circuitChanged();
  }

  function startWave() {
    if (G.phase !== "plan") return;
    G.phase = "battle";
    G.emitT = EMIT;
    say(`Level ${lvl().key} is moving. Registers leave the base every 2 seconds.`);
    updateHUD();
  }

  function checkCleared() {
    if (G.phase === "battle" && !G.boss && G.zombies.length === 0) {
      G.score += 25;
      G.phase = "cleared";
      const from = G.level, lastLevel = from === LEVELS.length - 1;
      say(lastLevel ? "Adventure 2 complete. The 3-qubit boss lines up again for another round."
        : `Level ${lvl().key} cleared. Level ${LEVELS[G.level + 1].key} is lining up.`);
      setTimeout(() => {
        if (G.phase !== "cleared" || G.level !== from) return;
        if (lastLevel) previewLevel(); else loadLevel(G.level + 1, true);
      }, 1500);
    }
  }

  // ---------- Entities ----------
  function emit() {
    G.regs.push({ id: G.nextId++, x: 0.25, prevX: 0.25, amps: basisState(0), next: 0, measured: false,
      idx: null, sym: null, basis: null, lanes: new Array(ND + 1).fill(true), broken: false, bossDone: false, flash: 0 });
  }

  function passColumns(r) {
    while (r.next < COLS.length && r.x >= COLS[r.next].u) {
      const c = COLS[r.next];
      if (c.kind === "L") {
        r.amps = c.apply(r.amps);
        if (c.blockEnd && G.rounds > 1) {
          for (let k = 1; k < G.rounds; k++) r.amps = grRound(r.amps, G.oracle, G.diff);
          G.loopFlash = 0.35;
        }
      } else {
        r.amps = applyCell(r.amps, cellOf(c));
      }
      r.next++;
    }
  }

  function measureReg(r) {
    const m = measureData(r.amps, G.basis);
    r.amps = m.amps; r.idx = m.d; r.basis = G.basis; r.measured = true; r.flash = 0.3;
    r.lanes[ANC] = false;
    r.sym = dataQs().map(q => (G.basis === "Z" ? String(bitOf(m.d, q)) : (bitOf(m.d, q) ? "-" : "+")));
    G.gateFlash = 0.25;
    G.lastOutcome = `${G.basis} basis → ${dataQs().reverse().map(q => sym(r.sym[q])).join("")}`;
  }

  function spawnZombie(type, lane) {
    if (G.phase !== "battle") { say("Test zombies can only be deployed while a level is moving."); return; }
    G.zombies.push({ id: G.nextId++, type, lane: Math.min(lane, ND - 1), x: LANE_END, prevX: LANE_END, hits: 0,
      alive: true, flash: 0, phase: Math.random() * 6, wait: 0 });
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
      ? "Measurement gate destroyed. Registers now stay quantum and travel the whole lane to x = 15."
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
    G.loopFlash = Math.max(0, G.loopFlash - dt);
    G.gateFlash = Math.max(0, G.gateFlash - dt);
    G.houseFlash = Math.max(0, G.houseFlash - dt);
    for (const s of G.splats) { s.life -= dt; s.r += dt * 50; }
    G.splats = G.splats.filter(s => s.life > 0);
    if (G.phase !== "battle" || G.paused) return;

    G.emitT += dt;
    while (G.emitT >= EMIT) { G.emitT -= EMIT; emit(); }

    for (const r of G.regs) {
      r.prevX = r.x; r.x += QSPEED * dt; r.flash = Math.max(0, r.flash - dt);
      passColumns(r);
      if (!r.measured && gateAlive() && r.x >= MEASURE_X) measureReg(r);
      if (r.x > LANE_END + 0.5 || r.lanes.slice(0, ND).every(l => !l)) r.dead = true;
    }

    for (const z of G.zombies) {
      if (!z.alive) continue;
      if (z.wait > 0) { z.wait -= dt; continue; }
      z.prevX = z.x; z.x -= ZSPEED * dt; z.flash = Math.max(0, z.flash - dt);
      if (z.type !== "Q" && gateAlive() && z.prevX >= MEASURE_X && z.x < MEASURE_X) {
        z.alive = false; splat(MEASURE_X, z.lane, "#ffd23f");
        damageGate(1, "A classical zombie hit the measurement gate and was destroyed");
        continue;
      }
      if (z.x <= 0.3) { z.alive = false; damageBase(1, `A ${z.type === "Q" ? "quantum" : "classical"} zombie reached the base`); }
    }

    const B = G.boss;
    if (B) {
      B.prevX = B.x; B.x -= BOSS_SPEED * dt; B.flash = Math.max(0, B.flash - dt);
      if (!B.hitGate && gateAlive() && B.x <= MEASURE_X) { B.hitGate = true; damageGate(2, "The boss rammed the measurement gate"); }
      if (B.x <= 0.3) { G.boss = null; damageBase(3, "The boss reached the base"); }
    }

    for (const r of G.regs) {
      if (!r.measured) continue;
      for (let q = 0; q < ND; q++) {
        if (!r.lanes[q]) continue;
        for (const z of G.zombies) {
          if (!z.alive || z.type === "Q" || z.lane !== q || z.wait > 0) continue;
          const cx = (r.prevX + r.x + z.prevX + z.x) / 4;
          if (cx < MEASURE_X - 0.05 || !crossed(r.prevX, r.x, z.prevX, z.x, RAD)) continue;
          r.lanes[q] = false; r.broken = true;
          if (r.sym[q] === z.type) {
            z.hits++; z.flash = 0.2; splat(z.x, q, "#b6f27a");
            if (z.hits >= 2) defeat(z, 10, `Two matching ${sym(z.type)} hits defeated a classical zombie on q${q}.`);
            else say(`Matching ${sym(z.type)} on q${q}. One more hit needed.`);
          } else {
            splat(z.x, q, "#d9d9d9");
            say(`Outcome ${sym(r.sym[q])} on q${q} did not match zombie ${sym(z.type)}.`);
          }
          break;
        }
      }
    }

    for (const r of G.regs) {
      if (r.measured) continue;
      const end = gateAlive() ? MEASURE_X : LANE_END;
      for (let q = 0; q < ND; q++) {
        if (!r.lanes[q]) continue;
        for (const z of G.zombies) {
          if (!z.alive || z.type !== "Q" || z.lane !== q || z.wait > 0) continue;
          const cx = (r.prevX + r.x + z.prevX + z.x) / 4;
          if (cx < COLS[0].u || cx > end || !crossed(r.prevX, r.x, z.prevX, z.x, RAD)) continue;
          const lab = blochLabel(r.amps, q);
          const m = measureQubitX(r.amps, q);
          r.amps = m.amps; r.lanes[q] = false; r.broken = true;
          if (m.plus) defeat(z, 20, `Qubit q${q} (${lab}) was found in |+⟩. Quantum zombie dissolved.`);
          else { splat(z.x, q, "#c9a7ff"); say(`Qubit q${q} (${lab}) was found in |−⟩. The quantum zombie shrugged it off.`); }
          break;
        }
      }
    }

    if (G.boss) {
      const Bs = G.boss;
      for (const r of G.regs) {
        if (!r.measured || r.bossDone || r.broken) continue;
        if (!crossed(r.prevX, r.x, Bs.prevX, Bs.x, RAD)) continue;
        r.bossDone = true;
        r.lanes = r.lanes.map(() => false);
        Bs.tests++;
        if (r.basis === "Z" && Bs.patterns.includes(r.idx)) {
          Bs.hits++; Bs.flash = 0.2;
          for (let q = 0; q < ND; q++) if (bitOf(r.idx, q)) splat(Bs.x, q, "#b6f27a");
          say(`Full-pattern hit ${ket(r.idx)}. Hits ${Bs.hits} of ${Bs.tests}.`);
        } else {
          say(`${r.basis === "Z" ? ket(r.idx) : "An X-basis outcome"} missed the boss. Hits ${Bs.hits} of ${Bs.tests}.`);
        }
        if (wilsonLow(Bs.hits, Bs.tests) > lvl().threshold) {
          G.score += 50;
          for (let q = 0; q < ND; q++) splat(Bs.x, q, "#ffd23f");
          say(`Boss defeated. ${Bs.hits} hits in ${Bs.tests} shots prove your hit rate beats the ${pct(lvl().threshold)} kill bar.`);
          G.boss = null;
          break;
        }
      }
    }

    G.regs = G.regs.filter(r => !r.dead);
    G.zombies = G.zombies.filter(z => z.alive);
    checkCleared();
  }

  // ---------- Duel ----------
  function duel() {
    if (G.phase !== "battle" || !G.boss) return;
    if (!gateAlive()) { say("The measurement gate is destroyed, so no duel is possible."); return; }
    let a = fullCircuit(G.oracle, G.diff, G.rounds);
    if (G.basis === "X") a = gHdata(a);
    const p = dataProbs(a), counts = new Array(DDIM).fill(0);
    for (let s = 0; s < 1000; s++) counts[sampleIndex(p)]++;
    const B = G.boss, need = Math.round(lvl().threshold * 1000);
    const hits = G.basis === "Z" ? B.patterns.reduce((s, i) => s + counts[i], 0) : 0;
    renderBars($("duelBars"), counts.map(c => c / 1000), counts.map(String), G.basis === "Z" ? B.patterns : []);
    if (hits > need) {
      G.score += 50;
      $("duelInfo").textContent = `${hits} of 1000 shots matched ${bossKet(B)}. More than ${need}: boss defeated.`;
      for (let q = 0; q < ND; q++) splat(B.x, q, "#ffd23f");
      G.boss = null;
      say("Duel won. The boss collapses.");
    } else {
      B.x -= 1.5;
      $("duelInfo").textContent = G.basis === "X"
        ? "X-basis outcomes are + and −, so none can match a lane pattern. The boss surges forward."
        : `Only ${hits} of 1000 shots matched. Needed more than ${need}. The boss surges forward.`;
      say("Duel lost. The boss surges forward.");
    }
    updateHUD();
  }

  // ---------- Circuit editing ----------
  function editCell(cell, lane) {
    const t = G.tool;
    if (t === "ERASE") { if (cell.ctrl) cell.ctrl = null; else if (lane < ND) cell.x[lane] = false; return true; }
    if (t === "X") {
      if (lane === ANC) { say("The ancilla's preparation is locked. Use ⊕ to target it."); return false; }
      if (cell.ctrl) cell.ctrl = null;
      cell.x[lane] = !cell.x[lane];
      return true;
    }
    const type = t;
    const qs = G.ctrlQs.slice();
    const same = cell.ctrl && cell.ctrl.type === type && cell.ctrl.qs.join() === qs.join();
    cell.x = new Array(ND).fill(false);
    cell.ctrl = same ? null : { type, qs };
    return true;
  }

  function gateName(ctrl) {
    const n = ctrl.qs.length;
    if (ctrl.type === "MCX") return n === 1 ? "CNOT" : n === 2 ? "Toffoli" : "C3X";
    if (ctrl.type === "MCSX") return n === 1 ? "C-√X" : n === 2 ? "CC-√X" : "C3-√X";
    if (ctrl.type === "MCS") return n === 1 ? "S" : n === 2 ? "CS" : "CCS";
    return n === 1 ? "Z" : n === 2 ? "CZ" : "CCZ";
  }

  function qiskitCode() {
    const all = dataQs().join(", ");
    const hAll = dataQs().map(q => `qc.h(${q});`).join(" ");
    const cells = [...G.oracle, ...G.diff];
    const usesSX = cells.some(c => c.ctrl && c.ctrl.type === "MCSX" && c.ctrl.qs.length > 1);
    const usesS = cells.some(c => c.ctrl && c.ctrl.type === "MCS" && c.ctrl.qs.length > 2);
    const cellLine = c => {
      if (c.ctrl) {
        const q = c.ctrl.qs, n = q.length;
        if (c.ctrl.type === "MCX") return n === 1 ? `qc.cx(${q[0]}, ${ANC});` : n === 2 ? `qc.ccx(${q[0]}, ${q[1]}, ${ANC});` : `qc.mcx([${q.join(", ")}], ${ANC});`;
        if (c.ctrl.type === "MCSX") return n === 1 ? `qc.csx(${q[0]}, ${ANC});` : `qc.append(new SXGate().control(${n}), [${q.join(", ")}, ${ANC}]);   // ${gateName(c.ctrl)}`;
        if (c.ctrl.type === "MCS") return n === 1 ? `qc.s(${q[0]});` : n === 2 ? `qc.cs(${q[0]}, ${q[1]});` : `qc.append(new SGate().control(2), [${q.join(", ")}]);   // CCS`;
        return n === 1 ? `qc.z(${q[0]});` : n === 2 ? `qc.cz(${q[0]}, ${q[1]});` : `qc.h(${q[2]}); qc.ccx(${q[0]}, ${q[1]}, ${q[2]}); qc.h(${q[2]});   // CCZ`;
      }
      return c.x.map((on, q) => (on ? `qc.x(${q});` : "")).filter(Boolean).join(" ") || null;
    };
    const imports = ["QuantumCircuit", ...(usesSX ? ["SXGate"] : []), ...(usesS ? ["SGate"] : [])].join(", ");
    const L = [`import { ${imports} } from 'qiskit';`, "",
      `// qubits ${all}: data lanes, qubit ${ANC}: ancilla j`,
      `const qc = QuantumCircuit(${ND + 1}, ${ND});`,
      `${hAll}   // init (locked)`,
      `qc.x(${ANC}); qc.h(${ANC});   // ancilla to |−⟩ (locked)`];
    for (let r = 0; r < G.rounds; r++) {
      L.push("", `// Grover round ${r + 1}: oracle`);
      G.oracle.forEach(c => { const s = cellLine(c); if (s) L.push(s); });
      L.push(`${hAll}   // locked`, "// diffuser middle");
      G.diff.forEach(c => { const s = cellLine(c); if (s) L.push(s); });
      L.push(`${hAll}   // locked`);
    }
    L.push("", `qc.h(${ANC});   // uncompute ancilla (locked)`,
      `qc.measure([${all}], [${all}]);${G.basis === "X" ? "   // add H on the data qubits first for X basis" : ""}`,
      "// next step: run on the pre-created Qollab backend, 1000 shots; keys read high qubit first.");
    return L.join("\n");
  }

  function circuitChanged() {
    $("qiskitCode").textContent = qiskitCode();
    let a = fullCircuit(G.oracle, G.diff, G.rounds);
    if (G.basis === "X") a = gHdata(a);
    const p = dataProbs(a);
    const pats = G.boss && G.basis === "Z" ? G.boss.patterns : [];
    $("predBars").classList.toggle("hidden", !G.showPred);
    renderBars($("predBars"), p, p.map(v => `${Math.round(v * 100)}%`), pats);
    if (!G.showPred) $("predInfo").textContent = "Prediction hidden. Turn it on to check your circuit, or find out in battle.";
    else if (G.boss) {
      const pb = pats.reduce((s, i) => s + p[i], 0), t = lvl().threshold;
      $("predInfo").textContent = G.basis === "X"
        ? "X basis: outcomes are + and −, so no shot can match the boss."
        : `P(boss) = ${(pb * 100).toFixed(1)}% with ${G.rounds} round${G.rounds > 1 ? "s" : ""}. ${pb > t + 1e-9 ? `Above the ${pct(t)} kill bar: your stream can prove it, or win a duel.` : `At or below the ${pct(t)} kill bar: this circuit cannot beat the boss.`}`;
    } else $("predInfo").textContent = "Output distribution of the planted circuit.";
    updateHUD();
  }

  // ---------- Panels ----------
  function buildBars(el, signed) {
    el.innerHTML = "";
    el.style.gridTemplateColumns = `repeat(${DDIM}, 1fr)`;
    for (let i = 0; i < DDIM; i++) {
      const bar = document.createElement("div");
      bar.className = "bar";
      bar.innerHTML = `<div class="track${signed ? " signed" : ""}"><div class="fill"></div><span class="val"></span></div><span class="lbl">${label(i)}</span>`;
      el.appendChild(bar);
    }
  }
  function renderBars(el, values, texts, highlight) {
    [...el.children].forEach((bar, i) => {
      bar.querySelector(".fill").style.height = `${Math.min(1, values[i]) * 100}%`;
      bar.querySelector(".val").textContent = texts[i];
      bar.classList.toggle("target", highlight.includes(i));
    });
  }
  function dataView(a) {
    const AB = 1 << ANC;
    const frames = [["|0⟩", [1, 0]], ["|1⟩", [0, 1]], ["|+⟩", [SQ, SQ]], ["|−⟩", [SQ, -SQ]]];
    let best = null;
    for (const [name, f] of frames) {
      const comp = Array.from({ length: DDIM }, (_, d) => cadd([f[0] * a[d][0], f[0] * a[d][1]], [f[1] * a[d | AB][0], f[1] * a[d | AB][1]]));
      const w = comp.reduce((s, v) => s + cabs2(v), 0);
      if (!best || w > best.w) best = { name, comp, w };
    }
    // Remove the global phase so the largest amplitude reads as a positive real number.
    const big = best.comp.reduce((m, v) => (cabs2(v) > cabs2(m) ? v : m), C(0));
    const n = Math.sqrt(cabs2(big)) || 1, g = [big[0] / n, -big[1] / n];
    best.comp = best.comp.map(v => cmul(v, g));
    return best;
  }
  function fmtAmp(v) {
    const [re, im] = v, e = 0.005;
    if (Math.abs(re) < e && Math.abs(im) < e) return "0";
    if (Math.abs(im) < e) return (re > 0 ? "+" : "−") + Math.abs(re).toFixed(2);
    if (Math.abs(re) < e) return (im > 0 ? "+" : "−") + Math.abs(im).toFixed(2) + "i";
    return `${re.toFixed(2)}${im < 0 ? "−" : "+"}${Math.abs(im).toFixed(2)}i`;
  }
  function blochLabel(a, q) {
    const { rz, rx, ry, r } = bloch(a, q);
    if (r < 0.98) return "entangled";
    if (rz > 0.98) return "|0⟩";
    if (rz < -0.98) return "|1⟩";
    if (rx > 0.98) return "|+⟩";
    if (rx < -0.98) return "|−⟩";
    if (ry > 0.98) return "|+i⟩";
    if (ry < -0.98) return "|−i⟩";
    return "ψ";
  }
  function updateMonitor() {
    const r = [...G.regs].reverse().find(x => x.next > 0);
    const a = r ? r.amps : prep(basisState(0));
    const v = dataView(a);
    [...$("ampBars").children].forEach((bar, i) => {
      const val = v.comp[i], fill = bar.querySelector(".fill");
      const real = Math.abs(val[1]) < 0.005, up = real ? val[0] >= 0 : true;
      fill.style.height = `${Math.sqrt(cabs2(val)) * 48}%`;
      fill.classList.toggle("neg", real && val[0] < -1e-9);
      fill.classList.toggle("cpx", !real);
      fill.style.bottom = up ? "50%" : "auto";
      fill.style.top = up ? "auto" : "50%";
      bar.querySelector(".val").textContent = fmtAmp(val);
    });
    const frame = `Data amplitudes with the ancilla in ${v.name} (${Math.round(v.w * 100)}% of the state).`;
    $("regInfo").textContent = !r ? `No register yet; showing the state after the locked prep. ${frame}`
      : r.measured ? `Register #${r.id} measured: ${G.lastOutcome}. ${frame}`
      : `Register #${r.id} at x = ${r.x.toFixed(2)}, ${r.next} of ${COLS.length} columns applied. ${frame}`;
    $("laneStates").style.gridTemplateColumns = `repeat(${ND + 1}, 1fr)`;
    $("laneStates").innerHTML = [...dataQs(), ANC].map(q => {
      const s = !r ? blochLabel(a, q) : !r.lanes[q] ? (q === ANC && r.measured ? "discarded" : "spent")
        : r.measured ? sym(r.sym[q]) : blochLabel(a, q);
      return `<div class="${q === ANC ? "anc" : ""}"><span>${q === ANC ? "j (ancilla)" : `q${q}`}</span><strong>${s}</strong></div>`;
    }).join("");
  }
  function syncControls() {
    document.querySelectorAll("#ctrlSel button").forEach(b => {
      const q = +b.dataset.q;
      b.classList.toggle("hidden", q >= ND);
      b.classList.toggle("on", G.ctrlQs.includes(q));
    });
    document.querySelectorAll("#laneSel button").forEach(b => {
      b.classList.toggle("hidden", +b.dataset.l >= ND);
      b.classList.toggle("on", +b.dataset.l === G.lane);
    });
    document.querySelectorAll("#rounds button").forEach(b => b.classList.toggle("on", +b.dataset.r === G.rounds));
    document.querySelectorAll("#levelSel button").forEach((b, i) => b.classList.toggle("on", i === G.level));
    document.querySelectorAll("#tools button").forEach(b => {
      b.classList.toggle("hidden", !lvl().tools.includes(b.dataset.tool));
      b.classList.toggle("on", b.dataset.tool === G.tool);
    });
    $("ctrlRow").classList.toggle("dim", G.tool === "X" || G.tool === "ERASE");
    $("textbook").classList.toggle("hidden", !lvl().textbook);
    $("toolHint").textContent = TOOL_HINT[G.tool];
    $("ctrlHint").textContent = `${ND} data qubits: pick 1 to ${ND} controls`;
  }
  function updateHUD() {
    $("levelChip").textContent = `${lvl().key} · ${G.level + 1} of ${LEVELS.length}`;
    $("baseHP").textContent = `${G.baseHP} / 10`;
    $("gateHP").textContent = gateAlive() ? `${G.gateHP} / 10` : "destroyed";
    $("score").textContent = G.score;
    $("killChip").textContent = `> ${Math.round(lvl().threshold * 1000)} / 1000`;
    $("bossChip").textContent = G.boss ? `${bossKet(G.boss)} · hits ${G.boss.hits}/${G.boss.tests}` : "none";
    $("pause").textContent = G.paused ? "Resume" : "Pause";
    $("basis").textContent = `Measurement basis: ${G.basis} (${G.basis === "Z" ? "0 / 1" : "+ / −"}), click to switch`;
    $("ready").classList.toggle("hidden", G.phase !== "plan");
    $("ready").textContent = `Ready: start level ${lvl().key}`;
    $("duel").disabled = G.phase !== "battle" || !G.boss || !gateAlive();
    document.querySelectorAll("[data-z]").forEach(b => { b.disabled = G.phase !== "battle"; });
    $("entityCount").textContent = `Registers ${G.regs.length} · Zombies ${G.zombies.length}`;
  }

  // ---------- Drawing ----------
  function rr(x, y, w, h, r) { ctx.beginPath(); ctx.roundRect(x, y, w, h, r); }
  function circle(x, y, r) { ctx.beginPath(); ctx.arc(x, y, r, 0, Math.PI * 2); }
  function diamond(x, y, r) { ctx.beginPath(); ctx.moveTo(x, y - r); ctx.lineTo(x + r, y); ctx.lineTo(x, y + r); ctx.lineTo(x - r, y); ctx.closePath(); }
  function text(s, x, y, size, color, align = "center") {
    ctx.font = `${size}px ${FONT}`; ctx.fillStyle = color; ctx.textAlign = align; ctx.textBaseline = "middle";
    ctx.fillText(s, x, y);
  }

  function drawLawn(t) {
    ctx.fillStyle = "#3d8a2a"; ctx.fillRect(0, 0, W, TOP);
    ctx.fillStyle = "#4c9c33";
    for (let x = 0; x < W; x += 34) { circle(x + 17, TOP - 2, 18); ctx.fill(); }
    for (let q = 0; q < ND; q++) {
      for (let u = 0; u < LANE_END; u++) {
        ctx.fillStyle = (q + u) % 2 ? "#6cba43" : "#7cc952";
        ctx.fillRect(xPos(u), TOP + q * LANE_H, xPos(u + 1) - xPos(u) + 1, LANE_H);
      }
    }
    ctx.fillStyle = "rgba(255,255,240,.12)";
    ctx.fillRect(xPos(0.4), TOP, xPos(MEASURE_X) - xPos(0.4), DATA_H);
    // Divider between the data lanes and the ancilla: a stone border with a dashed rule.
    const dy = TOP + DATA_H;
    ctx.fillStyle = "#9b8f7a"; ctx.fillRect(0, dy, W, GAP);
    ctx.fillStyle = "#bdb19a";
    for (let x = 0; x < W; x += 28) { rr(x + 2, dy + 2, 24, GAP - 4, 3); ctx.fill(); }
    ctx.save(); ctx.setLineDash([10, 6]); ctx.strokeStyle = "#f2b33d"; ctx.lineWidth = 3;
    ctx.beginPath(); ctx.moveTo(0, dy + GAP / 2); ctx.lineTo(W, dy + GAP / 2); ctx.stroke(); ctx.restore();
    // Ancilla pool lane
    const g = ctx.createLinearGradient(0, POOL_TOP, 0, BOTTOM);
    g.addColorStop(0, "#4fb6dd"); g.addColorStop(1, "#21709e");
    ctx.fillStyle = g; ctx.fillRect(0, POOL_TOP, W, POOL_H);
    ctx.strokeStyle = "rgba(255,255,255,.35)"; ctx.lineWidth = 2;
    for (let x = 110; x < W; x += 70) {
      ctx.beginPath(); ctx.arc(x + Math.sin(t + x) * 6, POOL_TOP + 20 + (x % 3) * 12, 10, Math.PI * 1.1, Math.PI * 1.9); ctx.stroke();
    }
    ctx.fillStyle = "#b98a4e"; ctx.fillRect(0, BOTTOM, W, H - BOTTOM);
    for (const u of [0, 5, 10, 15]) text(`x = ${u}`, xPos(u), H - 12, 12, "#5b3a1a");
    text("ancilla pool: no zombies, it only carries the phase kickback", xPos(9.5), laneY(ANC), 13, "rgba(255,255,255,.88)");
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
    text(`Base ${G.baseHP}`, 46, laneY(ANC), 15, "#fff");
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
  // Ancilla source: a lily pad with a golden diamond bud.
  function drawLily(x, y, t) {
    ctx.save(); ctx.translate(x, y + Math.sin(t * 1.5) * 1.5);
    ctx.beginPath(); ctx.ellipse(0, 6, 22, 9, 0, 0.35, Math.PI * 2 - 0.35); ctx.lineTo(0, 6); ctx.closePath();
    ctx.fillStyle = "#3f9e4a"; ctx.fill(); ctx.strokeStyle = "#1f5e27"; ctx.lineWidth = 2; ctx.stroke();
    diamond(0, -6, 11); ctx.fillStyle = "#ffcf4a"; ctx.fill(); ctx.strokeStyle = "#9a6a10"; ctx.stroke();
    ctx.restore();
  }

  function tile(cx, y, s, anc, lock) {
    if (anc) {
      diamond(cx, y, 17); ctx.fillStyle = "#ffe9a6"; ctx.fill(); ctx.lineWidth = 2.5; ctx.strokeStyle = "#9a6a10"; ctx.stroke();
      text(s, cx, y + 1, 14, "#7a4a00");
    } else {
      rr(cx - 14, y - 14, 28, 28, 6); ctx.fillStyle = "#dff6c8"; ctx.fill();
      ctx.lineWidth = 2.5; ctx.strokeStyle = "#2f7d1f"; ctx.stroke();
      text(s, cx, y + 1, 15, "#2f7d1f");
    }
    if (lock) { rr(cx + 6, y - 19, 9, 8, 2); ctx.fillStyle = "#8a5a2b"; ctx.fill(); }
  }
  function xTile(cx, y) {
    rr(cx - 14, y - 14, 28, 28, 6); ctx.fillStyle = "#ffb347"; ctx.fill();
    ctx.lineWidth = 2.5; ctx.strokeStyle = "#9a5b10"; ctx.stroke();
    text("X", cx, y + 1, 15, "#5b3a1a");
  }

  function drawCircuit(t) {
    const x1 = xPos(0.3), x2 = xPos(MEASURE_X);
    for (const q of dataQs()) {
      ctx.strokeStyle = "rgba(60,40,15,.45)"; ctx.lineWidth = 2;
      ctx.beginPath(); ctx.moveTo(x1, laneY(q)); ctx.lineTo(x2, laneY(q)); ctx.stroke();
      drawSprout(xPos(0.17), laneY(q) - 4, t);
      text(`q${q}`, xPos(0.17) - 2, laneY(q) + 38, 12, "#2a4a12");
    }
    ctx.strokeStyle = "#ffe9a6"; ctx.lineWidth = 2;
    ctx.beginPath(); ctx.moveTo(x1, laneY(ANC)); ctx.lineTo(xPos(COLS[10].u) + 18, laneY(ANC)); ctx.stroke();
    drawLily(xPos(0.17), laneY(ANC), t);
    text("j", xPos(0.17) - 26, laneY(ANC), 15, "#ffe9a6");

    COLS.forEach(c => {
      const cx = xPos(c.u);
      if (c.kind === "L") {
        for (const [q, s] of Object.entries(c.tiles)) tile(cx, laneY(+q), s, +q === ANC, true);
        return;
      }
      const cell = cellOf(c);
      if (cell.ctrl) {
        const qs = cell.ctrl.qs, ty = cell.ctrl.type;
        if (ty === "MCX" || ty === "MCSX") {
          ctx.strokeStyle = "#3b2412"; ctx.lineWidth = 3;
          ctx.beginPath(); ctx.moveTo(cx, laneY(Math.min(...qs))); ctx.lineTo(cx, laneY(ANC)); ctx.stroke();
          if (ty === "MCX") {
            circle(cx, laneY(ANC), 13); ctx.fillStyle = "#ffe9a6"; ctx.fill(); ctx.strokeStyle = "#9a6a10"; ctx.stroke();
            ctx.beginPath(); ctx.moveTo(cx - 13, laneY(ANC)); ctx.lineTo(cx + 13, laneY(ANC));
            ctx.moveTo(cx, laneY(ANC) - 13); ctx.lineTo(cx, laneY(ANC) + 13); ctx.stroke();
          } else {
            rr(cx - 17, laneY(ANC) - 14, 34, 28, 6); ctx.fillStyle = "#ffe9a6"; ctx.fill();
            ctx.lineWidth = 2.5; ctx.strokeStyle = "#9a6a10"; ctx.stroke();
            text("√X", cx, laneY(ANC) + 1, 14, "#7a4a00");
          }
          qs.forEach(q => { circle(cx, laneY(q), 7); ctx.fillStyle = "#3b2412"; ctx.fill(); });
        } else {
          const s = ty === "MCS" ? "S" : "Z", col = ty === "MCS" ? "#1f6f9a" : "#5b3a8c";
          if (qs.length === 1) {
            rr(cx - 14, laneY(qs[0]) - 14, 28, 28, 6); ctx.fillStyle = ty === "MCS" ? "#d2f1ff" : "#e9d7ff"; ctx.fill();
            ctx.lineWidth = 2.5; ctx.strokeStyle = col; ctx.stroke();
            text(s, cx, laneY(qs[0]) + 1, 15, col);
          } else {
            ctx.strokeStyle = col; ctx.lineWidth = 3;
            ctx.beginPath(); ctx.moveTo(cx, laneY(Math.min(...qs))); ctx.lineTo(cx, laneY(Math.max(...qs))); ctx.stroke();
            qs.forEach(q => { circle(cx, laneY(q), 7); ctx.fillStyle = col; ctx.fill(); });
            if (ty === "MCS") {
              const ym = (laneY(Math.min(...qs)) + laneY(Math.max(...qs))) / 2;
              rr(cx - 11, ym - 11, 22, 22, 5); ctx.fillStyle = "#d2f1ff"; ctx.fill(); ctx.lineWidth = 2; ctx.strokeStyle = col; ctx.stroke();
              text("S", cx, ym + 1, 13, col);
            }
          }
        }
      }
      if (c.reg === "D" && diffLocked()) {
        for (const q of dataQs()) { rr(cx + 6, laneY(q) - 21, 9, 8, 2); ctx.fillStyle = "#8a5a2b"; ctx.fill(); }
      }
      for (const q of dataQs()) {
        if (cell.ctrl && cell.ctrl.qs.includes(q)) continue;
        if (!cell.ctrl && cell.x[q]) { xTile(cx, laneY(q)); continue; }
        ctx.save(); ctx.setLineDash([4, 4]); ctx.lineWidth = 1.5; ctx.strokeStyle = "rgba(255,255,255,.7)";
        rr(cx - 11, laneY(q) - 11, 22, 22, 5); ctx.stroke(); ctx.restore();
      }
    });
    const bracket = (c1, c2, s) => {
      const a = xPos(COLS[c1].u) - 17, b = xPos(COLS[c2].u) + 17, y = 30;
      ctx.strokeStyle = "#c8402c"; ctx.lineWidth = 2;
      ctx.beginPath(); ctx.moveTo(a, y + 6); ctx.quadraticCurveTo(a, y, a + 8, y); ctx.lineTo(b - 8, y);
      ctx.quadraticCurveTo(b, y, b, y + 6); ctx.stroke();
      ctx.font = `13px ${FONT}`;
      const tw = ctx.measureText(s).width + 14;
      rr((a + b) / 2 - tw / 2, y - 20, tw, 18, 5); ctx.fillStyle = "rgba(255,248,220,.94)"; ctx.fill();
      text(s, (a + b) / 2, y - 11, 13, "#c8402c");
    };
    bracket(0, 1, "init");
    bracket(2, 9, G.loopFlash > 0 ? `looping ×${G.rounds}` : `Grover ×${G.rounds}`);
    bracket(10, 10, "undo");
    text("oracle", (xPos(COLS[2].u) + xPos(COLS[4].u)) / 2, BOTTOM + 16, 13, "#5b3a1a");
    text("diffuser", (xPos(COLS[6].u) + xPos(COLS[8].u)) / 2, BOTTOM + 16, 13, "#5b3a1a");
  }

  function drawGate(t) {
    const x = xPos(MEASURE_X);
    if (!gateAlive()) {
      ctx.fillStyle = "#7a5a35";
      for (const q of dataQs()) { ctx.save(); ctx.translate(x, laneY(q) + 30); ctx.rotate(0.5 - q * 0.4); ctx.fillRect(-12, -3, 24, 6); ctx.restore(); }
      text("gate destroyed", x, BOTTOM + 34, 12, "#8a2a1a");
      return;
    }
    ctx.save();
    if (G.gateFlash > 0) { ctx.shadowColor = "#fff3a0"; ctx.shadowBlur = 50 * G.gateFlash; }
    rr(x - 12, TOP - 2, 24, DATA_H + 4, 6); ctx.fillStyle = "#9a6a35"; ctx.fill();
    ctx.restore();
    ctx.lineWidth = 3; ctx.strokeStyle = "#5b3a1a"; rr(x - 12, TOP - 2, 24, DATA_H + 4, 6); ctx.stroke();
    for (const q of dataQs()) {
      const y = laneY(q);
      circle(x, y, 18); ctx.fillStyle = "#fff8dc"; ctx.fill(); ctx.strokeStyle = "#5b3a1a"; ctx.lineWidth = 2.5; ctx.stroke();
      const ang = -Math.PI / 2 + (G.gateFlash > 0 ? Math.sin(t * 40 + q) * 0.9 : 0.5 * Math.sin(t + q));
      ctx.beginPath(); ctx.moveTo(x, y + 5); ctx.lineTo(x + Math.cos(ang) * 12, y + 5 + Math.sin(ang) * 12);
      ctx.strokeStyle = "#c8402c"; ctx.lineWidth = 2.5; ctx.stroke();
    }
    rr(x - 52, BOTTOM + 6, 104, 22, 6); ctx.fillStyle = "#8a5a2b"; ctx.fill();
    text(`Measure ${G.basis} · ${G.gateHP} HP`, x, BOTTOM + 17, 12, "#fff4cf");
  }

  function drawRegs() {
    for (const r of G.regs) {
      const X = xPos(r.x);
      if (!r.measured) {
        const ent = [...dataQs(), ANC].filter(q => r.lanes[q] && bloch(r.amps, q).r < 0.98);
        if (ent.length >= 2) {
          ctx.save(); ctx.setLineDash([5, 4]); ctx.strokeStyle = "rgba(120,70,200,.85)"; ctx.lineWidth = 2.5;
          ctx.beginPath(); ctx.moveTo(X, laneY(Math.min(...ent))); ctx.lineTo(X, laneY(Math.max(...ent))); ctx.stroke(); ctx.restore();
        }
      }
      for (const q of [...dataQs(), ANC]) {
        if (!r.lanes[q]) continue;
        if (q === ANC) {
          ctx.save(); ctx.shadowColor = "#ffcf4a"; ctx.shadowBlur = 12;
          diamond(X, laneY(ANC), 9); ctx.fillStyle = "#ffcf4a"; ctx.fill();
          ctx.lineWidth = 2; ctx.strokeStyle = "#9a6a10"; ctx.stroke(); ctx.restore();
          continue;
        }
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
          if (r.x > COLS[10].u + 0.1) {
            const l = blochLabel(r.amps, q);
            text(l === "entangled" ? "ent" : l, X, y - 17, 12, "#3b1f6e");
          }
        }
      }
    }
  }

  function drawZombieBody(x, y, t, o) {
    const sw = o.still ? 0 : Math.sin(t * 2.6 + o.phase);
    ctx.save();
    ctx.globalAlpha = o.alpha;
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

  function drawZombies(t) {
    for (const z of G.zombies) {
      const x = xPos(z.x) - 8, y = laneY(z.lane) + 20, q = z.type === "Q";
      drawZombieBody(x, y, t, { scale: 0.82, alpha: 1, phase: z.phase, flash: z.flash > 0,
        still: G.phase === "plan" || z.wait > 0, coat: q ? "#6b4fa0" : "#6e5a43", skin: q ? "#c9b8ef" : "#9fc27a" });
      const card = q ? "Q |+⟩" : `${sym(z.type)}  ${z.hits}/2`;
      rr(x - 26, y - 90, 52, 18, 5); ctx.fillStyle = q ? "#5b3a8c" : "#8a5a2b"; ctx.fill();
      text(card, x, y - 81, 12, "#fff4cf");
    }
  }

  function drawBoss(t) {
    const B = G.boss;
    if (!B) return;
    // A body stands on every lane that is 1 in some pattern; lanes that differ between patterns are ghostly.
    const x = xPos(B.x);
    for (const q of dataQs()) {
      const bits = B.patterns.map(p => bitOf(p, q));
      if (!bits.some(Boolean)) continue;
      const ghost = !bits.every(Boolean);
      drawZombieBody(x, laneY(q) + 24, t, { scale: 1.05, alpha: ghost ? 0.5 : 1, phase: q,
        flash: B.flash > 0, still: G.phase === "plan", coat: "#4a3a5e", skin: "#a8d07f", crown: true });
    }
    // Card: the boss pattern, then a proof bar. The bar is the 3σ lower bound on your hit rate;
    // the boss falls when it passes the kill-bar mark. It never heals.
    const t0 = lvl().threshold, s = `Boss ${bossKet(B)}`, line = `hits ${B.hits} of ${B.tests} · prove more than ${pct(t0)}`;
    ctx.font = `14px ${FONT}`;
    const w1 = ctx.measureText(s).width;
    ctx.font = `12px ${FONT}`;
    const bw = Math.max(170, w1, ctx.measureText(line).width) + 24;
    const bx = Math.min(W - bw / 2 - 6, Math.max(bw / 2 + 6, x)), gx0 = bx - bw / 2 + 12, gw = bw - 24;
    rr(bx - bw / 2, 3, bw, 56, 7); ctx.fillStyle = "rgba(40,20,60,.92)"; ctx.fill();
    text(s, bx, 14, 14, "#efe6ff");
    rr(gx0, 26, gw, 8, 4); ctx.fillStyle = "#3a2a4a"; ctx.fill();
    rr(gx0, 26, gw * Math.max(0.02, wilsonLow(B.hits, B.tests)), 8, 4); ctx.fillStyle = "#b6f27a"; ctx.fill();
    ctx.fillStyle = "#ff6b5e"; ctx.fillRect(gx0 + gw * t0 - 1.5, 22, 3, 16);
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
    Q: { title: "Quantum zombie |+⟩", lines: [
      "It lives in the lab zone, before the gate.",
      "An unmeasured qubit in |+⟩ dissolves it;",
      "other states may bounce off."] },
    boss: { title: "Grover boss", lines: [
      "It hides on one secret bit pattern: a body stands on each 1.",
      "Mark that pattern with the oracle so measured bullets read it,",
      "and prove your hit rate beats the kill bar."] },
    pair: { title: "Superposed boss", lines: [
      "It hides in two patterns at once (the ghostly lane differs).",
      "Plain Grover gets stuck at exactly 50% here.",
      "Change the marking phase to break through."] }
  };
  const zombieKind = z => (z.type === "Q" ? "Q" : z.type === "+" || z.type === "-" ? "X" : "Z");
  const bossKind = () => (G.boss.patterns.length > 1 ? "pair" : "boss");
  const zombieAnchor = z => ({ x: xPos(z.x) - 8, top: laneY(z.lane) - 70, bottom: laneY(z.lane) + 46 });
  const bossSpan = () => ({ x: xPos(G.boss.x) - 52, y1: TOP + 8, y2: TOP + DATA_H - 8 });
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
    drawLawn(t); drawHouse(); drawCircuit(t); drawGate(t);
    drawZombies(t); drawBoss(t); drawRegs(); drawSplats(); drawCallouts(t);
    if (G.phase === "plan") {
      const a = xPos(6.0), b = xPos(11.4);
      rr(a, TOP + 8, b - a, 34, 8); ctx.fillStyle = "rgba(40,25,10,.8)"; ctx.fill();
      text(`Level ${lvl().key}: ${lvl().name}. Press Ready when your circuit is planted.`, (a + b) / 2, TOP + 25, 14, "#fff4cf");
    }
    if (G.paused || G.phase === "over") {
      ctx.fillStyle = "rgba(20,40,10,.55)"; ctx.fillRect(0, 0, W, H);
      text(G.phase === "over" ? "Base lost" : "Paused", W / 2, H / 2 - 10, 40, "#fff4cf");
      text(G.phase === "over" ? "Restart to play again" : "Resume when ready", W / 2, H / 2 + 26, 16, "#fff4cf");
    }
  }

  // ---------- Events ----------
  document.querySelectorAll("#tools button").forEach(b => b.addEventListener("click", () => {
    G.tool = b.dataset.tool;
    syncControls();
  }));
  document.querySelectorAll("#ctrlSel button").forEach(b => b.addEventListener("click", () => {
    const q = +b.dataset.q;
    const next = G.ctrlQs.includes(q) ? G.ctrlQs.filter(x => x !== q) : [...G.ctrlQs, q].sort();
    if (!next.length) { say("Keep at least one control qubit selected."); return; }
    G.ctrlQs = next;
    syncControls();
  }));
  document.querySelectorAll("#rounds button").forEach(b => b.addEventListener("click", () => {
    G.rounds = +b.dataset.r;
    syncControls();
    say(`The Grover block now runs ${G.rounds}×.`); circuitChanged();
  }));
  document.querySelectorAll("#laneSel button").forEach(b => b.addEventListener("click", () => {
    G.lane = +b.dataset.l; syncControls();
  }));
  document.querySelectorAll("#levelSel button").forEach((b, i) => b.addEventListener("click", () => {
    if (G.phase === "battle" && !confirm("Leave the level in progress?")) return;
    loadLevel(i, false);
  }));
  document.querySelectorAll("[data-z]").forEach(b => b.addEventListener("click", () => spawnZombie(b.dataset.z, G.lane)));

  $("textbook").addEventListener("click", () => {
    G.diff = textbookDiff();
    say(`Textbook diffuser planted: X on all data qubits, ${ND === 2 ? "CZ" : "CCZ"}, X on all.`); circuitChanged();
  });
  $("clearCircuit").addEventListener("click", () => {
    G.oracle = [emptyCell(), emptyCell(), emptyCell()];
    if (!diffLocked()) G.diff = [emptyCell(), emptyCell(), emptyCell()];
    say("Circuit cleared. Only the locked gates remain."); circuitChanged();
  });
  $("basis").addEventListener("click", () => {
    if (G.phase === "over") return;
    G.basis = G.basis === "Z" ? "X" : "Z"; say(`Measurement basis switched to ${G.basis}.`); circuitChanged();
  });
  $("showPred").addEventListener("change", e => { G.showPred = e.target.checked; circuitChanged(); });
  $("ready").addEventListener("click", startWave);
  $("duel").addEventListener("click", duel);
  $("pause").addEventListener("click", () => { if (G.phase === "battle") { G.paused = !G.paused; updateHUD(); } });
  $("restart").addEventListener("click", () => loadLevel(0, false));

  canvas.addEventListener("click", ev => {
    const rc = canvas.getBoundingClientRect();
    const x = (ev.clientX - rc.left) * (W / rc.width), y = (ev.clientY - rc.top) * (H / rc.height);
    let lane;
    if (y >= TOP && y < TOP + DATA_H) lane = Math.floor((y - TOP) / LANE_H);
    else if (y >= POOL_TOP && y <= BOTTOM) lane = ANC;
    else return;
    const c = COLS.find(cc => Math.abs(xPos(cc.u) - x) < 19);
    if (!c) return;
    if (c.kind === "L" || (c.reg === "D" && diffLocked())) { say("This gate is locked in place."); return; }
    const cell = cellOf(c);
    if (editCell(cell, lane)) {
      if (cell.ctrl) say(`${gateName(cell.ctrl)} planted on ${cell.ctrl.qs.map(q => `q${q}`).join(", ")}${cell.ctrl.type === "MCX" || cell.ctrl.type === "MCSX" ? " → j" : ""}.`);
      circuitChanged();
    }
  });

  // ---------- Tutorial storybook ----------
  // Pages shown before each level: a kid-friendly explanation, a "Try this" box and a collapsible
  // Physicist's corner. gt() draws a gate the way it looks on the lawn.
  const gt = (s, k) => `<span class="gt gt-${k}">${s}</span>`;
  const TUT = {
    intro: [{
      title: "The hide-and-seek zombie",
      kid: `<p>A giant boss zombie is hiding inside a <b>secret pattern</b> of 0s and 1s across your lanes. With 2 lanes there are 4 patterns: 00, 01, 10 and 11. To beat it you must fire that exact pattern at it, again and again.</p>
        <p>Guessing would be slow. Quantum computers have a famous trick for searching called <b>Grover's algorithm</b>. It makes the right answer louder and all the wrong answers quieter, at the same time.</p>`,
      try: ["The boss's pattern is written on its card. Your job is to build the machine that finds it.", "You can reopen these pages at any time with <b>How to play</b>."],
      phys: `<p>Grover's search finds one marked item among N using about (π/4)√N oracle calls, where a classical search needs about N/2 guesses on average: a quadratic speed-up. Here N = 4 or 8, small enough to watch every amplitude.</p>`
    }],
    levels: [
      [{
        title: "Step 1: everyone gets a raffle ticket",
        kid: `<p>The locked ${gt("H", "h")} gates at the start put the qubits into a <b>superposition of every pattern</b>, all with the same size. It's like giving each pattern one raffle ticket.</p>
          <p>If you measured right now, each of the 4 patterns would win 1 time in 4.</p>`,
        phys: `<p>H⊗H|00⟩ = |s⟩ = ½(|00⟩ + |01⟩ + |10⟩ + |11⟩): the uniform superposition with amplitude 1/√N on every basis state.</p>`
      }, {
        title: "Step 2: the oracle paints the boss",
        kid: `<p>The <b>oracle</b> secretly paints the boss's pattern with a <b>minus sign</b>. You can't see a sign by measuring, but the next step uses it.</p>
          <p>The painter is the helper qubit in the <b>pool</b>: the gold diamond in lane j. It sits in the state |−⟩, and when ${gt("⊕ onto j", "c")} flips it, the minus sign jumps back onto the pattern that caused the flip. This is called <b>phase kickback</b>.</p>`,
        try: [`${gt("⊕ onto j", "c")} only fires when <b>all</b> its control lanes are 1.`,
          `So: ${gt("X", "x")} on every lane where the boss pattern has a 0, then ${gt("⊕ onto j", "c")} with all controls, then the same ${gt("X", "x")} again to undo.`,
          `Example: boss |01⟩ means q1 = 0 and q0 = 1. Plant ${gt("X", "x")} on q1, then ${gt("⊕ onto j", "c")}, then ${gt("X", "x")} on q1.`],
        phys: `<p>X|−⟩ = −|−⟩, so a multi-controlled X onto the ancilla maps |x⟩|−⟩ → (−1)^f(x) |x⟩|−⟩ with f(x) = 1 only on the target. On the data register this is the oracle O = I − 2|w⟩⟨w|.</p>`
      }, {
        title: "Step 3: the diffuser makes it louder",
        kid: `<p>The <b>diffuser</b> flips every amplitude around the <b>average</b>. Thanks to its minus sign, the boss's pattern sits far below the average, so it bounces far above everyone else.</p>
          <p>With 4 patterns, a single round makes the boss's pattern come up <b>100%</b> of the time!</p>`,
        phys: `<p>D = 2|s⟩⟨s| − I is inversion about the mean. G = DO is a rotation by 2θ in the plane spanned by |w⟩ and |s⟩, where sin θ = √(M/N). For N = 4 and M = 1, θ = 30°, and one round reaches 3θ = 90°: exactly |w⟩.</p>`
      }, {
        title: "The duel: more than 500 out of 1000",
        kid: `<p>In this level the diffuser is already planted and <b>locked</b> for you. You only build the oracle.</p>
          <p>Press <b>Duel</b> to fire 1000 shots. If more than <b>500</b> land on the boss's pattern, it falls. Every register that reaches it also counts as a shot. The green bar on its card shows how sure the game is that you beat the kill bar; when it passes the red mark, the boss falls. It never heals.</p>`,
        try: ["Tick <b>Show predicted probabilities</b> to check your circuit before the battle.", "Press <b>Ready</b> to start the wave, then <b>Duel</b>."],
        phys: `<p>The boss never heals. Each register that reaches it is one shot: hit or miss. It falls when the 3σ (Wilson) lower bound on your hit rate passes the kill bar, or when a 1000-shot duel beats the bar. A perfect circuit (every shot a hit) proves more than 50% after 10 hits.</p>`
      }],
      [{
        title: "Build the diffuser yourself",
        kid: `<p>Now you plant the diffuser too. The recipe goes between the locked H walls:</p>
          <p>${gt("X", "x")} on every lane, then ${gt("Controlled-Z", "z")} on all lanes, then ${gt("X", "x")} on every lane again.</p>
          <p>The kill bar is higher: more than <b>700</b> out of 1000. Proving that takes lots of shots, so use the <b>Duel</b> to finish the boss.</p>`,
        phys: `<p>X⊗ⁿ · CZ · X⊗ⁿ puts a phase −1 on |0…0⟩ only. Sandwiched between Hadamard layers it equals −(2|s⟩⟨s| − I), the diffuser up to a global phase.</p>`
      }, {
        title: "More rounds isn't always better",
        kid: `<p>The <b>Repeat</b> buttons run the oracle and diffuser again. More sounds better, right? Not quite!</p>
          <p>Grover is like pushing a swing. Push at the right moment and it flies high. Keep pushing and it comes back down. With 4 patterns: 1 round gives <b>100%</b>, 2 rounds <b>25%</b>, 3 rounds <b>25%</b>, 4 rounds <b>100%</b> again.</p>`,
        try: ["Try each Repeat value with predictions on and watch the bars jump."],
        phys: `<p>After k rounds P = sin²((2k+1)θ). With θ = 30° the angle runs 90°, 150°, 210°, 270°: Grover is a rotation, so overshooting is real.</p>`
      }],
      [{
        title: "A boss in two places at once",
        kid: `<p>This boss is <b>superposed</b>: it hides in <b>two</b> patterns at the same time, like (|01⟩ + |11⟩)/√2. A shot on either pattern hits it.</p>
          <p>Try the ordinary Grover recipe and something strange happens: the chance is stuck at exactly <b>50%</b>, no matter how many rounds you run! The kill bar is more than 500, so 50% isn't enough.</p>`,
        phys: `<p>With M = 2 marked states out of N = 4, sin θ = 1/√2 so θ = 45°, and sin²((2k+1)·45°) = 1/2 for every k. Standard Grover just rotates back and forth across the same point.</p>`
      }, {
        title: "The quarter-turn trick",
        kid: `<p>Ordinary Grover paints a minus sign, which is a <b>half turn</b>. Here we paint a <b>quarter turn</b>, called <b>i</b>, instead.</p>
          <p>${gt("√X onto j", "c")} does it in the oracle: the pool helper turns its half-flip into an i. ${gt("Controlled-S", "r")} does it in the diffuser. With both, <b>one round finds the boss 100%</b> of the time!</p>`,
        try: [`The two boss patterns share one lane's bit. Use that lane as the only control of ${gt("√X onto j", "c")}, with ${gt("X", "x")} around it if that bit is 0.`,
          `Diffuser: ${gt("X", "x")} on all, ${gt("Controlled-S", "r")} on q0 and q1, ${gt("X", "x")} on all.`,
          "Using the trick in only one of the two places leaves you at 50%: both are needed."],
        phys: `<p>√X|−⟩ = i|−⟩, so the oracle multiplies marked states by i, and X·CS·X between Hadamards gives I − (1 − i)|s⟩⟨s|. For M/N = 1/2 one iteration cancels the unmarked amplitude exactly. This is phase-matched (exact) Grover search: G. L. Long, Phys. Rev. A 64, 022307 (2001).</p>`
      }],
      [{
        title: "A bigger garden: 8 patterns",
        kid: `<p>Three lanes means <b>8 patterns</b>, so the boss is harder to find. The recipe is the same, just bigger: the oracle needs ${gt("⊕ onto j", "c")} with all three controls, and the diffuser ${gt("Controlled-Z", "z")} on all three lanes.</p>`,
        try: ["Use the same X-sandwich trick on the lanes where the boss has a 0.", "<b>Plant textbook diffuser</b> is available as a shortcut."],
        phys: `<p>sin θ = 1/√8, so θ ≈ 20.7°. The 3-control gates are a C3X onto the ancilla and a CCZ on the data qubits.</p>`
      }, {
        title: "Count your rounds",
        kid: `<p>The kill bar here is more than <b>900</b> out of 1000. One round gives about <b>78%</b>: not enough! Two rounds give about <b>94.5%</b>: that wins. Three rounds overshoot to 33%, and four fall to about 1%.</p>`,
        try: ["Set Repeat to 2 before you duel."],
        phys: `<p>The best round count is about (π/4)√N − ½ ≈ 1.7, so 2 rounds. The success curve sin²((2k+1)θ) peaks and falls again because Grover is a rotation, not a steady climb.</p>`
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
  // Each level's pages open once per session, the first time it loads (with the welcome page before level a).
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
    last = t; G.t = t;
    update(dt);
    draw(t);
    G.monT += dt;
    if (G.monT > 0.15) { G.monT = 0; updateMonitor(); updateHUD(); }
    requestAnimationFrame(frame);
  }

  window.QPVZ = { G, fullCircuit, dataProbs, duel, startWave, loadLevel, update, circuitChanged };
  fitCanvas();
  window.addEventListener("resize", fitCanvas);
  loadLevel(0, false);
  updateMonitor();
  requestAnimationFrame(frame);
})();
