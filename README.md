# Quantum Battlefield: a quantum twist on Plants vs. Zombies

**Your bullets are qubits. Your plants are quantum gates. The zombies obey quantum mechanics.**

▶ **Play in your browser:** https://Sriharsha-13.github.io/quantum-battlefield/

Every shot fired across the lawn is a live quantum state. It stays in superposition, carrying amplitudes and phases, until the wooden measurement gate collapses it into an ordinary classical pea. The gates you plant shape that state, and each zombie falls only to the state you actually built. Every outcome is sampled from the exact quantum state by the Born rule.

Built by Sriharsha for the **Qollab Hackathon 2026**.

## A Qollab Hackathon project

Quantum Battlefield was created for the **Qollab Hackathon 2026**. A big thank-you to the Qollab team for organising the hackathon and for the chance to build this game.

🧪 **Run it in the Qollab workspace:** [Quantum Battlefield on Qollab](https://qollab.xyz/p/01a123a4-1b9e-79ae-b26d-b35089b05541)

## Adventures

| Adventure | Levels | Physics you play with |
|---|---|---|
| **1. Qubit basics** | A. First qubit · B. Rotations · C. Entangled boss | Superposition, phase, measurement bases, rotations, Bell states, an entanglement witness |
| **2. Grover's search** | a. Mark the boss · b. Build the diffuser · c. Superposed boss · d. Three lanes | Oracles, phase kickback, the diffuser, over-rotation, phase-matched Grover for two targets |
| **3. Ground-state hunters** | a. Commuting ladder · b. Trotter–Suzuki | Variational energy estimation (VQE), Trotter–Suzuki steps, adiabatic state preparation, shot noise |

Each level opens with a short illustrated tutorial for kids, with a "Physicist's corner" for the full story. New zombies and bosses are introduced on the lawn the first time they appear, and every level shows its equivalent **Qiskit** circuit.

Levels unlock in order, and you earn up to 3 stars per level. Progress is saved in your browser. Judges and teachers can tick **Unlock all levels** on any adventure map.

For the physical intuition behind each adventure, and the references, see [project_page.md](project_page.md).

## Run it locally

No install or build step. Download or clone the repository and open `index.html` in a browser.

```
git clone https://github.com/Sriharsha-13/quantum-battlefield.git
```

## How the files are organised

```
index.html                 launcher: adventure tiles and level maps
launcher/                  launcher styles and scripts, progress and stars
1st adventure/             code.js · index.html · style.css · play.html
2nd adventure/             code.js · index.html · style.css · play.html
3rd adventure/             code.js · index.html · style.css · play.html
build-play-pages.ps1       regenerates each play.html
project_page.md            the Qollab project write-up
```

Each adventure folder keeps the three files the Qollab "JS / Qiskit" editor expects: `code.js`, an `index.html` body fragment and `style.css`. Those files can be pasted into Qollab unchanged. `play.html` wraps them into a full web page for the launcher. After editing an adventure's `index.html`, regenerate the wrappers:

```
powershell -ExecutionPolicy Bypass -File build-play-pages.ps1
```
