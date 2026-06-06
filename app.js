// Mathematical Constant Constants
const SQRT_2PI = Math.sqrt(2 * Math.PI);

// Abramowitz & Stegun Approximation for Error Function
function erf(x) {
    const a1 =  0.254829592;
    const a2 = -0.284496736;
    const a3 =  1.421413741;
    const a4 = -1.453152027;
    const a5 =  1.061405429;
    const p  =  0.3275911;

    const sign = (x < 0) ? -1 : 1;
    x = Math.abs(x);

    const t = 1.0 / (1.0 + p * x);
    const y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);

    return sign * y;
}

// Cumulative Distribution Function for Standard Normal
function normalCDF(z) {
    return 0.5 * (1.0 + erf(z / Math.sqrt(2.0)));
}

// Lanczos approximation for ln(Gamma(z))
function logGamma(z) {
    const g = 7;
    const p = [
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916283859,
        12.507381424447072,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7
    ];
    if (z < 0.5) {
        return Math.log(Math.PI) - Math.log(Math.sin(Math.PI * z)) - logGamma(1.0 - z);
    }
    z -= 1.0;
    let x = p[0];
    for (let i = 1; i < g + 2; i++) {
        x += p[i] / (z + i);
    }
    let t = z + g + 0.5;
    return 0.5 * Math.log(2.0 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x);
}

// Natural logarithm of n!
function logFactorial(n) {
    if (n <= 1) return 0;
    return logGamma(n + 1);
}

// Binomial PMF: P(H = k) = (N choose k) * p^k * (1-p)^(N-k)
function binomialPMF(k, N, p) {
    if (k < 0 || k > N) return 0;
    if (p === 0) return k === 0 ? 1 : 0;
    if (p === 1) return k === N ? 1 : 0;
    
    const lnP = logFactorial(N) - logFactorial(k) - logFactorial(N - k) 
              + k * Math.log(p) + (N - k) * Math.log(1 - p);
    return Math.exp(lnP);
}

// Binomial Tail Probability: P(H > X) = sum_{k=X+1}^{N} P(H=k)
// Optimised in O(min(tail_size, N)) with early termination
function binomialTailProbability(X, N, p) {
    const mean = N * p;
    let sum = 0;
    if (X >= mean) {
        // Sum right tail (X+1 to N), terms decrease
        for (let k = X + 1; k <= N; k++) {
            const pmf = binomialPMF(k, N, p);
            if (pmf < 1e-30) break;
            sum += pmf;
        }
        return sum;
    } else {
        // Sum left tail (0 to X), terms decrease as k -> 0
        for (let k = X; k >= 0; k--) {
            const pmf = binomialPMF(k, N, p);
            if (pmf < 1e-30) break;
            sum += pmf;
        }
        return 1.0 - sum;
    }
}

// Box-Muller transform for standard Normal distribution sampling
function randomNormal(mean, stdDev) {
    let u = 0, v = 0;
    while(u === 0) u = Math.random();
    while(v === 0) v = Math.random();
    let num = Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * v);
    return num * stdDev + mean;
}

// Fast approximate Binomial sampler using Normal distribution
function sampleBinomial(N, p) {
    const mean = N * p;
    const stdDev = Math.sqrt(N * p * (1 - p));
    let val = Math.round(randomNormal(mean, stdDev));
    return Math.max(0, Math.min(N, val));
}

// --- App State & UI Elements ---
let simActive = false;
let simInterval = null;
let chartInstance = null;

// Simulation parameters
let N = 10000;
let K = 10000;
let X = 5200;
let bias = 0.5;
let method = 'monte-carlo';

// Simulation counters
let totalTrialsRun = 0;
let successTrialsRun = 0;
let simulatedOutcomes = []; // Array of length N+1

// UI Elements
const numFlipsRange = document.getElementById('num-flips-range');
const numFlipsInput = document.getElementById('num-flips');
const numTrialsRange = document.getElementById('num-trials-range');
const numTrialsInput = document.getElementById('num-trials');
const thresholdRange = document.getElementById('threshold-x-range');
const thresholdInput = document.getElementById('threshold-x');
const threshValDisplay = document.getElementById('thresh-val-display');

const runBtn = document.getElementById('run-btn');
const cancelBtn = document.getElementById('cancel-btn');
const progressContainer = document.getElementById('sim-progress-container');
const progressBar = document.getElementById('sim-progress-bar');
const progressText = document.getElementById('sim-progress-text');

// Stats Displays
const cltProbDisplay = document.getElementById('clt-prob');
const binomProbDisplay = document.getElementById('binom-prob');
const empProbDisplay = document.getElementById('emp-prob');
const empRatioDisplay = document.getElementById('emp-ratio');
const zScoreDisplay = document.getElementById('z-score');

// Math Display Toggles
const tabButtons = document.querySelectorAll('.tab-btn');
const tabPanes = document.querySelectorAll('.tab-pane');

// Chart Toggles
const toggleNormal = document.getElementById('toggle-normal-curve');
const toggleBinomial = document.getElementById('toggle-binom-curve');

// Math equation details element bindings
const mathN = document.getElementById('math-n');
const mathP = document.getElementById('math-p');
const mathQ = document.getElementById('math-q');
const mathMeanCalc = document.getElementById('math-mean-calc');
const mathMean = document.getElementById('math-mean');
const mathVarCalc = document.getElementById('math-var-calc');
const mathVar = document.getElementById('math-var');
const mathSd = document.getElementById('math-sd');
const mathX = document.getElementById('math-x');
const mathXcc = document.getElementById('math-x-cc');
const mathXraw = document.getElementById('math-x-raw');
const mathMeanRaw1 = document.getElementById('math-mean-raw1');
const mathSdRaw1 = document.getElementById('math-sd-raw1');
const mathZraw = document.getElementById('math-z-raw');
const mathXccRaw = document.getElementById('math-x-cc-raw');
const mathMeanRaw2 = document.getElementById('math-mean-raw2');
const mathSdRaw2 = document.getElementById('math-sd-raw2');
const mathZcc = document.getElementById('math-z-cc');
const mathZrawVal = document.getElementById('math-z-raw-val');
const mathZrawVal2 = document.getElementById('math-z-raw-val2');
const mathProbRaw = document.getElementById('math-prob-raw');
const mathZccVal = document.getElementById('math-z-cc-val');
const mathZccVal2 = document.getElementById('math-z-cc-val2');
const mathProbCc = document.getElementById('math-prob-cc');
const mathXbinom = document.getElementById('math-x-binom');
const mathBinomResult = document.getElementById('math-binom-result');
const mathCcResultCompare = document.getElementById('math-cc-result-compare');

// Set up UI linkages (Range slider synced to numeric input)
function syncInputs(range, number, updateCallback) {
    range.addEventListener('input', () => {
        number.value = range.value;
        updateCallback();
    });
    number.addEventListener('change', () => {
        // Clamp values
        let val = parseInt(number.value);
        const min = parseInt(number.min);
        const max = parseInt(number.max);
        if (isNaN(val)) val = min;
        val = Math.max(min, Math.min(max, val));
        number.value = val;
        range.value = val;
        updateCallback();
    });
}

// Sync inputs
syncInputs(numFlipsRange, numFlipsInput, handleFlipsChange);
syncInputs(numTrialsRange, numTrialsInput, handleParamsChange);
syncInputs(thresholdRange, thresholdInput, handleParamsChange);

function handleFlipsChange() {
    // Threshold can't be larger than N
    N = parseInt(numFlipsInput.value);
    thresholdRange.max = N;
    thresholdInput.max = N;
    if (parseInt(thresholdInput.value) > N) {
        thresholdInput.value = N;
        thresholdRange.value = N;
    }
    handleParamsChange();
}

function handleParamsChange() {
    N = parseInt(numFlipsInput.value);
    K = parseInt(numTrialsInput.value);
    X = parseInt(thresholdInput.value);
    threshValDisplay.textContent = X;
    updateCalculations();
}

// Coin bias button configuration
document.querySelectorAll('.bias-btn').forEach(btn => {
    btn.addEventListener('click', (e) => {
        document.querySelectorAll('.bias-btn').forEach(b => b.classList.remove('active'));
        btn.classList.add('active');
        bias = parseFloat(btn.dataset.bias);
        updateCalculations();
    });
});

// Simulation method radio selection
document.querySelectorAll('.method-option').forEach(option => {
    option.addEventListener('click', () => {
        document.querySelectorAll('.method-option').forEach(opt => opt.classList.remove('active'));
        option.classList.add('active');
        const radio = option.querySelector('input[type="radio"]');
        radio.checked = true;
        method = radio.value;
    });
});

// Tab setup
tabButtons.forEach(btn => {
    btn.addEventListener('click', () => {
        tabButtons.forEach(b => b.classList.remove('active'));
        tabPanes.forEach(p => p.classList.remove('active'));
        
        btn.classList.add('active');
        const tabId = btn.dataset.tab;
        document.getElementById(tabId).classList.add('active');
    });
});

// Chart Toggles
toggleNormal.addEventListener('change', () => updateChartCurves());
toggleBinomial.addEventListener('change', () => updateChartCurves());

// Math Notation helper: round to readable scientific notation or decimal
function formatScientific(val) {
    if (val === 0) return "0.000%";
    if (val < 0.0001) {
        return val.toExponential(4);
    }
    return (val * 100).toFixed(5) + "%";
}

function formatRawDecimal(val) {
    if (val < 0.000001 && val > 0) {
        return val.toExponential(4);
    }
    return val.toFixed(7);
}

// Calculate statistical parameters and display formulas
function updateCalculations() {
    // 1. Core Params
    const p = bias;
    const q = 1 - p;
    const mean = N * p;
    const variance = N * p * q;
    const stdDev = Math.sqrt(variance);
    
    // Z-scores
    const zRaw = (X - mean) / stdDev;
    const zCc = ((X + 0.5) - mean) / stdDev;
    
    // Probabilities
    const probRawVal = 1.0 - normalCDF(zRaw);
    const probCcVal = 1.0 - normalCDF(zCc);
    const probBinomVal = binomialTailProbability(X, N, p);

    // Update main stats panel
    cltProbDisplay.textContent = formatScientific(probCcVal);
    binomProbDisplay.textContent = formatScientific(probBinomVal);
    zScoreDisplay.textContent = zCc.toFixed(4);

    // Update LaTeX values in math card
    mathN.textContent = N.toLocaleString();
    mathP.textContent = p;
    mathQ.textContent = q.toFixed(2);
    mathMeanCalc.textContent = `${N.toLocaleString()} \\cdot ${p}`;
    mathMean.textContent = mean.toLocaleString(undefined, {maximumFractionDigits: 2});
    mathVarCalc.textContent = `${N.toLocaleString()} \\cdot ${p} \\cdot ${q.toFixed(2)}`;
    mathVar.textContent = variance.toLocaleString(undefined, {maximumFractionDigits: 2});
    mathSd.textContent = stdDev.toFixed(4);

    mathX.textContent = X.toLocaleString();
    mathXcc.textContent = (X + 0.5).toLocaleString();
    mathXraw.textContent = X.toLocaleString();
    mathMeanRaw1.textContent = mean.toLocaleString(undefined, {maximumFractionDigits: 2});
    mathSdRaw1.textContent = stdDev.toFixed(4);
    mathZraw.textContent = zRaw.toFixed(4);

    mathXccRaw.textContent = (X + 0.5).toLocaleString();
    mathMeanRaw2.textContent = mean.toLocaleString(undefined, {maximumFractionDigits: 2});
    mathSdRaw2.textContent = stdDev.toFixed(4);
    mathZcc.textContent = zCc.toFixed(4);

    mathZrawVal.textContent = zRaw.toFixed(2);
    mathZrawVal2.textContent = zRaw.toFixed(2);
    mathProbRaw.textContent = formatRawDecimal(probRawVal);

    mathZccVal.textContent = zCc.toFixed(2);
    mathZccVal2.textContent = zCc.toFixed(2);
    mathProbCc.textContent = formatRawDecimal(probCcVal);

    mathXbinom.textContent = X.toLocaleString();
    mathBinomResult.textContent = formatRawDecimal(probBinomVal);
    mathCcResultCompare.textContent = formatRawDecimal(probCcVal);

    // Trigger MathJax re-render if loaded
    if (window.MathJax && typeof MathJax.typesetPromise === 'function') {
        MathJax.typesetPromise();
    }

    // Refresh charts
    initChartStructure();
}

// Chart bin config
function getChartBinConfig() {
    const mean = N * bias;
    const stdDev = Math.sqrt(N * bias * (1 - bias));
    
    // Choose dynamic bounds around the mean (usually 4.5 standard deviations)
    const bounds = 4.5;
    let xMin = Math.max(0, Math.floor(mean - bounds * stdDev));
    let xMax = Math.min(N, Math.ceil(mean + bounds * stdDev));
    
    // Determine optimal bin width to keep ~40-50 bins on chart
    const targetBins = 40;
    const range = xMax - xMin;
    let binWidth = Math.max(1, Math.round(range / targetBins));
    
    // Align bounds to binWidth multiples
    xMin = Math.floor(xMin / binWidth) * binWidth;
    xMax = Math.ceil(xMax / binWidth) * binWidth;
    
    const binsCount = Math.ceil((xMax - xMin) / binWidth) + 1;
    
    return { xMin, xMax, binWidth, binsCount };
}

// Initialize the chart template
function initChartStructure() {
    if (chartInstance) {
        chartInstance.destroy();
    }

    if (typeof Chart === 'undefined') {
        const wrapper = document.querySelector('.chart-wrapper');
        if (wrapper) {
            wrapper.innerHTML = `<div style="text-align: center; color: var(--accent-pink); padding: 2rem; border: 1px dashed var(--accent-pink); border-radius: var(--border-radius-md);">
                <strong>Chart.js not loaded</strong><br>
                <span style="font-size: 0.8rem; color: var(--text-secondary); display: block; margin-top: 0.5rem;">
                    Please ensure chart.js is present in the project folder and page is refreshed.
                </span>
            </div>`;
        }
        return;
    }

    const ctx = document.getElementById('distribution-chart').getContext('2d');
    const { xMin, xMax, binWidth, binsCount } = getChartBinConfig();
    
    // Construct labels and data templates
    const labels = [];
    const simulatedData = [];
    const normalData = [];
    const binomialData = [];
    const barBackgrounds = [];

    const mean = N * bias;
    const variance = N * bias * (1 - bias);
    const stdDev = Math.sqrt(variance);

    for (let i = 0; i < binsCount; i++) {
        const binStart = xMin + i * binWidth;
        const binEnd = binStart + binWidth;
        const binCenter = binStart + binWidth / 2;

        // Label representing the range
        if (binWidth === 1) {
            labels.push(binStart.toString());
        } else {
            labels.push(`${binStart}-${binEnd - 1}`);
        }

        // Initialize empty simulated frequency
        simulatedData.push(0);

        // Theoretical Normal Height for this bin (integrated over the bin width)
        // Prob(binStart <= Y < binEnd) = Phi((binEnd - mean)/stdDev) - Phi((binStart - mean)/stdDev)
        const zStart = (binStart - mean) / stdDev;
        const zEnd = (binEnd - mean) / stdDev;
        const normalProb = normalCDF(zEnd) - normalCDF(zStart);
        normalData.push(normalProb * K); // scale by total K trials

        // Theoretical Binomial Height
        let binomialProb = 0;
        // Sum binomial PMF inside this bin
        for (let k = binStart; k < binEnd; k++) {
            binomialProb += binomialPMF(k, N, bias);
        }
        binomialData.push(binomialProb * K);

        // Highlight bar colors if the bin is above the threshold X
        if (binStart > X) {
            barBackgrounds.push('rgba(255, 8, 68, 0.45)'); // Bright glowing neon pink
        } else if (binEnd - 1 > X) {
            barBackgrounds.push('rgba(255, 100, 100, 0.35)'); // Mixed boundary
        } else {
            barBackgrounds.push('rgba(0, 242, 254, 0.3)'); // Cool cyan
        }
    }

    chartInstance = new Chart(ctx, {
        type: 'bar',
        data: {
            labels: labels,
            datasets: [
                {
                    label: 'Simulated Outcomes',
                    data: simulatedData,
                    backgroundColor: barBackgrounds,
                    borderColor: barBackgrounds.map(c => c.replace('0.3', '0.6').replace('0.45', '0.8')),
                    borderWidth: 1.5,
                    barPercentage: 0.95,
                    categoryPercentage: 0.95,
                    order: 3
                },
                {
                    label: 'Theoretical Normal (CLT)',
                    data: normalData,
                    type: 'line',
                    borderColor: 'rgba(0, 242, 254, 1)',
                    backgroundColor: 'transparent',
                    borderWidth: 2.5,
                    pointRadius: 0,
                    tension: 0.35,
                    hidden: !toggleNormal.checked,
                    order: 1
                },
                {
                    label: 'Theoretical Binomial',
                    data: binomialData,
                    type: 'line',
                    borderColor: 'rgba(155, 81, 224, 1)',
                    borderDash: [5, 4],
                    backgroundColor: 'transparent',
                    borderWidth: 2.0,
                    pointRadius: 0,
                    tension: 0.35,
                    hidden: !toggleBinomial.checked,
                    order: 2
                }
            ]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            plugins: {
                legend: {
                    labels: {
                        color: '#9aa0c0',
                        font: { family: 'Plus Jakarta Sans', size: 11 }
                    }
                },
                tooltip: {
                    backgroundColor: 'rgba(12, 12, 30, 0.95)',
                    titleColor: '#fff',
                    bodyColor: '#b3b3d9',
                    borderColor: 'rgba(255, 255, 255, 0.08)',
                    borderWidth: 1,
                    titleFont: { family: 'Outfit', weight: 'bold' },
                    bodyFont: { family: 'Plus Jakarta Sans' },
                    callbacks: {
                        label: function(context) {
                            const val = context.parsed.y;
                            const pct = ((val / K) * 100).toFixed(4) + "%";
                            return `${context.dataset.label}: ${Math.round(val).toLocaleString()} trials (${pct})`;
                        }
                    }
                }
            },
            scales: {
                x: {
                    grid: { color: 'rgba(255, 255, 255, 0.03)' },
                    ticks: {
                        color: '#676c8c',
                        font: { family: 'Fira Code', size: 9 },
                        maxTicksLimit: 15
                    },
                    title: {
                        display: true,
                        text: 'Number of Heads',
                        color: '#9aa0c0',
                        font: { family: 'Plus Jakarta Sans', size: 12, weight: 'bold' }
                    }
                },
                y: {
                    grid: { color: 'rgba(255, 255, 255, 0.03)' },
                    ticks: {
                        color: '#676c8c',
                        font: { family: 'Fira Code', size: 9 }
                    },
                    title: {
                        display: true,
                        text: 'Frequency (Trials count)',
                        color: '#9aa0c0',
                        font: { family: 'Plus Jakarta Sans', size: 12, weight: 'bold' }
                    }
                }
            }
        }
    });
}

function updateChartCurves() {
    if (!chartInstance) return;
    chartInstance.data.datasets[1].hidden = !toggleNormal.checked;
    chartInstance.data.datasets[2].hidden = !toggleBinomial.checked;
    chartInstance.update('none'); // silent update
}

// Update the chart with current simulated outcome counts
function updateChartData() {
    if (!chartInstance) return;

    const { xMin, binWidth, binsCount } = getChartBinConfig();
    const binnedCounts = new Array(binsCount).fill(0);

    // Aggregate simulated outcomes into bins
    for (let h = 0; h < simulatedOutcomes.length; h++) {
        const count = simulatedOutcomes[h];
        if (count > 0) {
            // Find appropriate bin
            let binIndex = Math.floor((h - xMin) / binWidth);
            if (binIndex >= 0 && binIndex < binsCount) {
                binnedCounts[binIndex] += count;
            }
        }
    }

    // Update chart simulated dataset
    chartInstance.data.datasets[0].data = binnedCounts;
    
    // Scale theoretical curves to matches current actual total run
    const currentTotal = Math.max(1, totalTrialsRun);
    const mean = N * bias;
    const stdDev = Math.sqrt(N * bias * (1 - bias));

    const normalData = [];
    const binomialData = [];

    for (let i = 0; i < binsCount; i++) {
        const binStart = xMin + i * binWidth;
        const binEnd = binStart + binWidth;

        // Scale by current total run
        const zStart = (binStart - mean) / stdDev;
        const zEnd = (binEnd - mean) / stdDev;
        const normalProb = normalCDF(zEnd) - normalCDF(zStart);
        normalData.push(normalProb * currentTotal);

        let binomialProb = 0;
        for (let k = binStart; k < binEnd; k++) {
            binomialProb += binomialPMF(k, N, bias);
        }
        binomialData.push(binomialProb * currentTotal);
    }

    chartInstance.data.datasets[1].data = normalData;
    chartInstance.data.datasets[2].data = binomialData;
    
    chartInstance.update('none'); // update without animations for speed
}

// Reset the simulation metrics
function resetSimulationState() {
    totalTrialsRun = 0;
    successTrialsRun = 0;
    simulatedOutcomes = new Array(N + 1).fill(0);
    
    empProbDisplay.textContent = "0.00000%";
    empRatioDisplay.textContent = `0 / ${K.toLocaleString()} trials`;
    empProbDisplay.classList.remove('text-glow');

    initChartStructure();
}

// Run the simulation
function startSimulation() {
    if (simActive) return;
    
    simActive = true;
    runBtn.disabled = true;
    cancelBtn.disabled = false;
    
    // Disable parameters during simulation run
    numFlipsRange.disabled = true;
    numFlipsInput.disabled = true;
    numTrialsRange.disabled = true;
    numTrialsInput.disabled = true;
    thresholdRange.disabled = true;
    thresholdInput.disabled = true;
    document.querySelectorAll('.bias-btn').forEach(btn => btn.disabled = true);
    
    progressContainer.classList.remove('hidden');
    progressBar.style.width = '0%';
    progressText.textContent = 'Initializing simulation...';

    resetSimulationState();

    if (method === 'binomial-sample') {
        runFastSimulation();
    } else {
        runMonteCarloSimulation();
    }
}

// Fast simulation using Box-Muller normal approximation
function runFastSimulation() {
    const chunkCount = 10;
    const chunkSize = Math.ceil(K / chunkCount);
    let chunkIndex = 0;

    function processChunk() {
        if (!simActive) return;

        const start = chunkIndex * chunkSize;
        const end = Math.min(K, start + chunkSize);

        for (let i = start; i < end; i++) {
            const heads = sampleBinomial(N, bias);
            simulatedOutcomes[heads]++;
            totalTrialsRun++;
            if (heads > X) {
                successTrialsRun++;
            }
        }

        // Update progress UI
        const pct = (totalTrialsRun / K) * 100;
        progressBar.style.width = `${pct}%`;
        progressText.textContent = `Generating samples: ${totalTrialsRun.toLocaleString()} / ${K.toLocaleString()}`;

        // Update statistics
        const empiricalProb = successTrialsRun / totalTrialsRun;
        empProbDisplay.textContent = (empiricalProb * 100).toFixed(5) + "%";
        empRatioDisplay.textContent = `${successTrialsRun.toLocaleString()} / ${totalTrialsRun.toLocaleString()} trials`;
        if (successTrialsRun > 0) empProbDisplay.classList.add('text-glow');

        // Refresh visualization
        updateChartData();

        chunkIndex++;
        if (chunkIndex < chunkCount && totalTrialsRun < K) {
            simInterval = setTimeout(processChunk, 10);
        } else {
            endSimulation();
        }
    }

    processChunk();
}

// Detailed Monte Carlo Simulation (flipping virtual coins)
function runMonteCarloSimulation() {
    // Dynamic batch size to ensure stable 60fps (~1M operations per frame)
    const opsPerFrame = 1200000;
    const batchSize = Math.max(1, Math.floor(opsPerFrame / N));
    
    function processBatch() {
        if (!simActive) return;

        const trialsLeft = K - totalTrialsRun;
        const currentBatchSize = Math.min(batchSize, trialsLeft);

        for (let t = 0; t < currentBatchSize; t++) {
            let heads = 0;
            // Actually simulate each individual coin flip
            for (let f = 0; f < N; f++) {
                if (Math.random() < bias) {
                    heads++;
                }
            }

            simulatedOutcomes[heads]++;
            totalTrialsRun++;
            if (heads > X) {
                successTrialsRun++;
            }
        }

        // Update progress UI
        const pct = (totalTrialsRun / K) * 100;
        progressBar.style.width = `${pct}%`;
        progressText.textContent = `Flipping coins: ${totalTrialsRun.toLocaleString()} / ${K.toLocaleString()}`;

        // Update statistics
        const empiricalProb = successTrialsRun / totalTrialsRun;
        empProbDisplay.textContent = (empiricalProb * 100).toFixed(5) + "%";
        empRatioDisplay.textContent = `${successTrialsRun.toLocaleString()} / ${totalTrialsRun.toLocaleString()} trials`;
        if (successTrialsRun > 0) empProbDisplay.classList.add('text-glow');

        // Refresh visualization
        updateChartData();

        if (totalTrialsRun < K) {
            simInterval = requestAnimationFrame(processBatch);
        } else {
            endSimulation();
        }
    }

    simInterval = requestAnimationFrame(processBatch);
}

function stopSimulation() {
    if (!simActive) return;
    
    if (method === 'binomial-sample') {
        clearTimeout(simInterval);
    } else {
        cancelAnimationFrame(simInterval);
    }
    
    progressText.textContent = `Cancelled at ${totalTrialsRun.toLocaleString()} trials.`;
    endSimulation(true);
}

function endSimulation(cancelled = false) {
    simActive = false;
    runBtn.disabled = false;
    cancelBtn.disabled = true;

    // Enable parameters
    numFlipsRange.disabled = false;
    numFlipsInput.disabled = false;
    numTrialsRange.disabled = false;
    numTrialsInput.disabled = false;
    thresholdRange.disabled = false;
    thresholdInput.disabled = false;
    document.querySelectorAll('.bias-btn').forEach(btn => btn.disabled = false);

    if (!cancelled) {
        progressText.textContent = `Completed! ${totalTrialsRun.toLocaleString()} trials simulated.`;
    }
}

// Event bindings
runBtn.addEventListener('click', startSimulation);
cancelBtn.addEventListener('click', stopSimulation);

// Load MathJax for typesetting
const mathjaxScript = document.createElement('script');
mathjaxScript.src = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js";
mathjaxScript.async = true;
document.head.appendChild(mathjaxScript);

// Start initial setup on page load
try {
    handleParamsChange();
} catch (error) {
    const errorDiv = document.createElement('div');
    errorDiv.style.position = 'fixed';
    errorDiv.style.top = '0';
    errorDiv.style.left = '0';
    errorDiv.style.width = '100%';
    errorDiv.style.background = '#ff0844';
    errorDiv.style.color = '#fff';
    errorDiv.style.padding = '1rem';
    errorDiv.style.zIndex = '99999';
    errorDiv.style.fontFamily = 'monospace';
    errorDiv.style.fontSize = '0.9rem';
    errorDiv.style.boxShadow = '0 4px 20px rgba(0,0,0,0.5)';
    errorDiv.innerHTML = '<strong>Initialization Error:</strong> ' + error.message + ' <br> <span style="font-size:0.75rem;">Stack: ' + error.stack + '</span>';
    document.body.insertBefore(errorDiv, document.body.firstChild);
}
