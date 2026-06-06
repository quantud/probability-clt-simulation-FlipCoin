# Probability & CLT Simulation — Flip Coin

An interactive web application demonstrating the **Central Limit Theorem (CLT)** through coin flip simulations.

## Features

- **Monte Carlo Simulation** — simulates individual coin flips in real time
- **Binomial Sampler** — fast statistical sampling using the Box-Muller transform
- **Analytical Solutions** — CLT normal approximation (with/without continuity correction) vs. exact binomial probability
- **Interactive Chart** — live histogram of outcomes with overlaid theoretical Normal and Binomial curves
- **Adjustable Parameters** — flips per trial (N), number of trials (K), heads threshold (X), coin bias (p)
- **Math Walkthrough** — step-by-step derivation of the Z-score and probability calculations using LaTeX rendered by MathJax

## Files

| File | Description |
|---|---|
| `index.html` | Main page — layout, inline script, Chart.js reference |
| `style.css` | Glassmorphism dark-mode design system |
| `app.js` | Standalone JS (math functions + simulation engine) |
| `chart.js` | Bundled Chart.js v4 (local, no CDN dependency) |
| `verify_math.ps1` | PowerShell script to verify analytical formulas independently |

## Usage

Open `index.html` directly in a browser — no build step or server required.

## Math

For N coin flips with bias p:
- μ = N·p
- σ = √(N·p·q)
- Z = (X − μ) / σ
- P(Heads > X) ≈ 1 − Φ(Z)  (CLT approximation)
- P(Heads > X) = Σ C(N,k) p^k (1−p)^(N−k)  (exact binomial)

## Tech Stack

- Vanilla HTML, CSS, JavaScript
- Chart.js v4 (bundled locally)
- MathJax 3 (CDN, for LaTeX rendering)
- Google Fonts: Outfit, Plus Jakarta Sans, Fira Code

---
© 2026 Coin Flip CLT Project
