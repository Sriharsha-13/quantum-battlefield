// Quantum Battlefield launcher glue for an adventure page. The adventure's own code.js is
// untouched: this opens the level named in ?level=N (0-based), locks the in-game level buttons
// the player hasn't reached yet, and awards stars when the game reports a level as cleared.
(function () {
  const P = window.QBProgress;
  const adv = +document.body.dataset.adv;
  const api = window.QPVZ1 || window.QPVZ || window.QPVZ3;
  const G = api && api.G;
  const buttons = Array.from(document.querySelectorAll("#levelSel button"));

  function syncLocks() {
    buttons.forEach((b, i) => {
      const open = P.unlocked(adv, i);
      b.disabled = !open;
      b.classList.toggle("qb-locked", !open);
      b.title = open ? "" : "Clear the level before this one to unlock it";
    });
  }
  function syncUrl(i) {
    const url = new URL(location.href);
    if (url.searchParams.get("level") === String(i)) return;
    url.searchParams.set("level", i);
    history.replaceState(null, "", url);
  }

  syncLocks();
  let want = parseInt(new URLSearchParams(location.search).get("level"), 10);
  if (!(want >= 0 && want < buttons.length)) want = 0;
  while (want > 0 && !P.unlocked(adv, want)) want--;
  if (want > 0) buttons[want].click();
  syncUrl(want);

  if (!G) return;
  let award = -1;
  setInterval(() => {
    if (G.phase === "cleared") {
      if (award !== G.level) {
        award = G.level;
        P.record(adv, G.level, P.starsFor(G.baseHP));
        syncLocks();
      }
    } else {
      award = -1;
    }
    syncUrl(G.level);
  }, 250);
})();
