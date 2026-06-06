function Erf($x) {
    $a1 = 0.254829592
    $a2 = -0.284496736
    $a3 = 1.421413741
    $a4 = -1.453152027
    $a5 = 1.061405429
    $p = 0.3275911

    $sign = 1
    if ($x -lt 0) { $sign = -1 }
    $x_abs = [Math]::Abs($x)

    $t = 1.0 / (1.0 + $p * $x_abs)
    $y = 1.0 - (((((($a5 * $t + $a4) * $t) + $a3) * $t + $a2) * $t + $a1) * $t * [Math]::Exp(-$x_abs * $x_abs))

    return $sign * $y
}

function NormalCDF($z) {
    return 0.5 * (1.0 + (Erf ($z / [Math]::Sqrt(2.0))))
}

function LogGamma($z) {
    $g = 7
    $p = @(
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916283859,
        12.507381424447072,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7
    )
    if ($z -lt 0.5) {
        return [Math]::Log([Math]::PI) - [Math]::Log([Math]::Sin([Math]::PI * $z)) - (LogGamma (1.0 - $z))
    }
    $z = $z - 1.0
    $x = $p[0]
    for ($i = 1; $i -lt ($g + 2); $i++) {
        $x += $p[$i] / ($z + $i)
    }
    $t = $z + $g + 0.5
    return 0.5 * [Math]::Log(2.0 * [Math]::PI) + ($z + 0.5) * [Math]::Log($t) - $t + [Math]::Log($x)
}

function LogFactorial($n) {
    if ($n -le 1) { return 0.0 }
    return LogGamma ($n + 1)
}

function BinomialPMF($k, $N, $p) {
    if ($k -lt 0 -or $k -gt $N) { return 0.0 }
    $lnP = (LogFactorial $N) - (LogFactorial $k) - (LogFactorial ($N - $k)) + $k * [Math]::Log($p) + ($N - $k) * [Math]::Log(1.0 - $p)
    return [Math]::Exp($lnP)
}

function BinomialTailProbability($X, $N, $p) {
    $mean = $N * $p
    $sum = 0.0
    if ($X -ge $mean) {
        for ($k = $X + 1; $k -le $N; $k++) {
            $pmf = BinomialPMF $k $N $p
            if ($pmf -lt 1e-30) { break }
            $sum += $pmf
        }
        return $sum
    } else {
        for ($k = $X; $k -ge 0; $k--) {
            $pmf = BinomialPMF $k $N $p
            if ($pmf -lt 1e-30) { break }
            $sum += $pmf
        }
        return 1.0 - $sum
    }
}

$N = 10000
$p = 0.5
$X = 5200

$mean = $N * $p
$variance = $N * $p * (1.0 - $p)
$stdDev = [Math]::Sqrt($variance)

$zRaw = ($X - $mean) / $stdDev
$pCltRaw = 1.0 - (NormalCDF $zRaw)

$zCc = (($X + 0.5) - $mean) / $stdDev
$pCltCc = 1.0 - (NormalCDF $zCc)

$pBinom = BinomialTailProbability $X $N $p

Write-Host "=== Verification of Analytical Solutions (PowerShell) ==="
Write-Host "Parameters: N = $N, p = $p, Threshold X = $X"
Write-Host "Mean (mu) = $mean"
Write-Host "Variance (sigma^2) = $variance"
Write-Host "Std Dev (sigma) = $stdDev"
Write-Host ""
Write-Host "--- Central Limit Theorem (Normal Approx) ---"
Write-Host "Without Continuity Correction:"
Write-Host "  Z-Score = $($zRaw.ToString('F4'))"
Write-Host "  P(H > 5200) = $($pCltRaw.ToString('E8')) ($(($pCltRaw * 100).ToString('F6'))%)"
Write-Host "With Continuity Correction:"
Write-Host "  Z-Score = $($zCc.ToString('F4'))"
Write-Host "  P(H > 5200) = $($pCltCc.ToString('E8')) ($(($pCltCc * 100).ToString('F6'))%)"
Write-Host ""
Write-Host "--- Exact Binomial Distribution ---"
Write-Host "  P(H > 5200) = $($pBinom.ToString('E8')) ($(($pBinom * 100).ToString('F6'))%)"
Write-Host ""
Write-Host "--- Observations ---"
$diff = [Math]::Abs($pCltCc - $pBinom)
Write-Host "Abs difference between CLT (with CC) and Exact Binomial: $($diff.ToString('E8'))"
Write-Host "Relative error: $( (($diff / $pBinom) * 100).ToString('F4') )%"
