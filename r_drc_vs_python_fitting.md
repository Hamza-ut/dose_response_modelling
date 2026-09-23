Reasoning log for why `fit_ll4` (in `drc/logistic.py`) exists, and why each
initial-guess choice was made the way it was.

R's `drc` package is the gold standard these are compared against — its results are
treated as ground truth throughout.

# Initial guesses

Initial guesses are the most important part of fitting a 4-parameter logistic (4PL)
curve with `scipy.optimize.curve_fit`. If they're bad, the whole fit can fail to
converge, or converge to a numerically meaningless answer.

## c (lower) & d (upper)

**R_DRC**: takes the observed `y_data` min/max and pads them slightly outward (0.1% of
the response range), so `c` sits just below the true minimum and `d` sits just above the
true maximum.

```python
y_min, y_max = y_data.min(), y_data.max()
pad = 0.001 * (y_max - y_min)
c_lower_initial = y_min - pad
d_upper_initial = y_max + pad
```

This uses every data point, which is good for normal, clean data. But it has a real
weakness: a single contaminated/outlier reading directly becomes the guess. Tested on
real data with one bad reading injected: `c` was dragged from a sane `~1.0` to `-20`,
essentially becoming the outlier itself.

**Hamza_DRC**: uses the 5th and 95th percentile of `y_data` instead of the raw min/max.

```python
c_lower_initial = np.percentile(y_data, 5)
d_upper_initial = np.percentile(y_data, 95)
```

This throws away the most extreme ~10% of points for this specific calculation, but is
far more resistant to a single bad reading — the same outlier test that dragged R's `c`
to `-20` only moved percentile's `c` to `-10.5`, and on clean data the two methods
produce nearly identical guesses anyway (`0.28` vs `0.22` on real ryegrass data). The
practical consequence is bigger than the guess itself: when both guesses were carried
all the way through to a final `curve_fit` result, R's outlier-corrupted guess caused a
complete numerical blowup (`c ≈ -7,900,000`, `e = 0.0`), while percentile's guess still
converged to a sane, if imperfect, fit (`b=18.1` — high enough to warrant scrutiny, but
finite and plausible).

**Verdict: percentile.** The excluded ~10% of points don't hurt in practice — they're
also the points nearest the asymptotes, which are the least informative for estimating
slope anyway, since the curve is nearly flat right at the asymptotes.

## b (slope) & e(ED50)

**R_DRC**: derives both `b` _and_ `e` directly from the data via a single linear
regression, instead of guessing either. This works because of the 4PL formula's own
algebra. Starting from

$$y = c + \frac{d-c}{1+(x/e)^b}$$

and rearranging step by step (subtract `c`, flip the fraction, subtract 1, take the
log of both sides) gives:

$$\log\!\left(\frac{d-y}{y-c}\right) = b\log(x) - b\log(e)$$

That's the equation of a straight line in `(log(x), Z)` space, where
`Z = log((d-y)/(y-c))` — with **slope = b** and **intercept = -b·log(e)**. So instead of
guessing `b`, compute `Z` for every real data point, fit a line through
`(log(x), Z)`, and read `b` straight off the slope — this isn't an arbitrary transform
that happens to look linear, it's the exact algebraic inverse of the model being fit.
If the data didn't actually follow a 4PL shape, this wouldn't produce a clean line at
all.

`e` falls out of the _same_ line for free, and there's a concrete reason the intercept
specifically is where it lives: **`e` is defined as the dose where the response sits
exactly halfway between `c` and `d`** (that's what "half-maximal effective dose" means).
Plug `x = e` into the original formula and `y` comes out to exactly `(c+d)/2` — so at
that one specific dose, `d - y` and `y - c` are equal, their ratio is exactly `1`, and
`log(1) = 0`. Plug `x = e` into the fitted line (`Z = b·log(x) + intercept`) and that
same `Z = 0` shows up on the left: `0 = b·log(e) + intercept`, which rearranges straight
to `intercept = -b·log(e)`, i.e. `e = exp(-intercept/slope)`. In other words, the
intercept isn't an abstract algebraic leftover — it's telling you exactly where the
fitted line crosses zero, and that crossing point _is_ the EC50 by definition. No
separate calculation, no separate look at the data, just reading a second meaningful
number off the line you already fit for `b`.

```python
positive_mask = x_data > 0        # log(x) is undefined at x <= 0
log_x = np.log(x_data[positive_mask])
transformed_y = np.log(
    (d_upper_initial - y_data[positive_mask])
    / (y_data[positive_mask] - c_lower_initial)
)
regression = linregress(log_x, transformed_y)
b_slope_initial = regression.slope
e_inflection_initial = np.exp(-regression.intercept / regression.slope)
```

Proof `b` actually tracks the real slope, not noise: fit synthetic data generated from
a known `b=3.5` — a fixed, data-blind guess of `1.0` is off by `2.5`; this regression
guess is off by only `0.38`. (Early versions of this project used a hardcoded
`b_slope_initial = 1.0` for every dataset, on the assumption that a 45° slope is a
reasonable generic default — abandoned for exactly this reason: it's the same number
regardless of what the data says, so it can't tell a gentle curve from a steep one.)

Proof `e` does the same: on real ryegrass data, `median(x_data)` (the old approach)
guessed `2.815`, `exp(-intercept/slope)` guessed `2.901` — the final converged fit
landed at `3.058`, so the regression guess was closer (`0.157` off vs `0.243` off). The
gap gets much bigger once the tested doses aren't evenly spread around the true EC50:
with a dose series covering `1` to `200` but weighted toward the low end (true EC50=50),
`median(x_data)` guessed `20.0` — more than half the true value — while the regression
guessed `39.88`, much closer. `curve_fit` recovered the correct final answer (`~51`)
from _both_ starting points in that case, but that's `curve_fit`'s own refinement doing
extra work to correct a poor starting guess — exactly the kind of slack that runs out
under the outlier scenario in the c/d section above, where a bad starting point made
the difference between a sane fit and a complete numerical blowup. (Old approach:
`e_ed50_initial = np.median(x_data)` — just "the middle dose tested," computed without
looking at the response data at all, so it can't tell whether the doses tested were
actually centered on the real transition or not.)

**Hamza_DRC**: adopts R's regression trick for both `b` and `e`, but feeds it `c`/`d`
from the percentile method above instead of R's raw min/max — same reasoning as the c/d
section: outlier-robustness without giving up the adaptiveness of a real, data-derived
estimate. Using the same regression for both parameters, rather than regression for `b`
and `median(x_data)` for `e`, also keeps the two guesses internally consistent — they
describe the same fitted line, not two unrelated readings of the data stitched together.

```python
c_lower_initial = np.percentile(y_data, 5)
d_upper_initial = np.percentile(y_data, 95)

positive_mask = (
    (y_data > c_lower_initial) & (y_data < d_upper_initial) & (x_data > 0)
)

log_x = np.log(x_data[positive_mask])
transformed_y = np.log(
    (d_upper_initial - y_data[positive_mask])
    / (y_data[positive_mask] - c_lower_initial)
)
regression = linregress(log_x, transformed_y)
b_slope_initial = regression.slope
e_inflection_initial = np.exp(-regression.intercept / regression.slope)
```

**The mask needs to check more than R's version does.** R's `positive_mask` only
excludes `x <= 0`, which is enough for R's min/max-based `c`/`d`, since every point is
guaranteed to sit inside `[c, d]` by construction there. Percentile-based `c`/`d` is
different: it _guarantees_ some points fall outside `[c_lower_initial,
d_upper_initial]` (that's what "5th/95th percentile" means) — for those points,
`(y - c_lower_initial)` or `(d_upper_initial - y)` is negative, so `log()` of the ratio
is `NaN`, and since `NaN` isn't skipped, it would poison the _entire_ regression, not
just those rows. So the mask has to check all three conditions.

**Added error checking.** After the `[c, d]` trimming, a very small dataset can be
left with too few points for a meaningful regression: `linregress` needs at least 2
points to run at all, but 2 points is a degenerate "fit" (zero residual degrees of
freedom — a line through 2 points always "fits perfectly," which says nothing about
reliability), and with only 1 point it silently returns `b = nan` instead of a real
error. The guard requires at least 4 points to survive (2 degrees of freedom: 4
points − 2 parameters being fit ≥ 2), raising a clear error naming exactly how many
points survived instead of continuing on with a meaningless or `NaN` result:

```python
n_used = positive_mask.sum()
if n_used < 4:
    raise ValueError(
        f"Only {n_used} point(s) remain after excluding non-positive doses "
        f"and points outside the 5th-95th percentile range; need at least "
        f"4 for a meaningful slope regression. Dataset has {len(y_data)} "
        f"points total."
    )
```

**Known limitation, shared by both R's and Hamza's version, not something either fixes:**
the regression estimate is only as good as the dose range actually collected. If the
tested doses don't extend far enough to approach the true asymptotes (common for gentle,
slow-plateauing curves), or don't have enough points near the actual inflection point
(common for very steep curves, where the transition happens over a narrow slice of the
dose range), the `b`/`e` estimates can still be systematically biased — regardless of
the c/d strategy used. This isn't an outlier problem, it's a "does the experiment's dose
range actually contain enough information" problem.

**Verdict:** adopt R's regression-based approach for both `b` and `e`, combined with
Hamza's percentile-based `c`/`d`.

## Optimization algorithm

Once the initial guesses are in hand, both R and scipy still have to actually search
for the minimum of the same objective (the sum of squared residuals between the curve
and the data). *How* they search is not the same, and it's worth being precise about
why, since it's tempting to assume R's choice is the "correct" one simply because a
biostatistician wrote it.

**R_DRC**: `drm()` fits by calling R's general-purpose optimizer `optim()`, and its
default `method` (confirmed from `drc`'s own `drmc()`/`drmOpt.R` source) is
**`"BFGS"`** — a quasi-Newton method that treats the objective as an arbitrary scalar
function to minimize. It only ever looks at the objective's value and gradient, and
builds up an approximate picture of the curvature from previous steps; it has no
built-in notion that the objective happens to be a sum of squares at all. This
generality is a deliberate, necessary tradeoff on R's side: `drm()` also fits
non-Gaussian response types (binomial, Poisson, etc.) through the same
objective-minimization machinery, so it needs one optimizer that works for *any*
log-likelihood, not one specialized for least-squares. If bounds are supplied, `drmOpt.R`
switches to `"L-BFGS-B"` instead — still the same general-purpose family, just with box
constraints bolted on.

**Hamza_DRC**: `scipy.optimize.curve_fit` picks between two solvers that are
*specialized* for least-squares problems specifically, and unlike R, the pick is
automatic based on whether `bounds=` is passed:

- No bounds (`use_bounds=False`, scipy's own default for an unconstrained call):
  **Levenberg-Marquardt** (`method='lm'`, via MINPACK's `leastsq`). It directly uses the
  residuals' Jacobian — how each individual data point's error changes as the
  parameters move — to take targeted steps toward the minimum, instead of inferring
  curvature indirectly the way BFGS does.
- With bounds (`use_bounds=True`): scipy switches to **Trust Region Reflective**
  (`method='trf'`) automatically, because MINPACK's `lm` has no way to respect a box
  constraint. `trf` is still a least-squares specialist (it uses the same Jacobian
  structure as `lm`), it just restricts each step to a region it currently trusts and
  reflects off the bounds instead of crossing them. It also computes a usable covariance
  matrix via a pseudoinverse when the Jacobian is rank-deficient (poorly-identified
  parameters), where `lm` just returns `inf`.

```python
fit_kwargs = {"p0": initial_guess, "maxfev": 5000}
if use_bounds:
    lower_bounds = [0, np.min(y_data) - 1, np.min(y_data), 0]
    upper_bounds = [100, np.max(y_data) + 1, np.max(y_data), np.max(x_data)]
    fit_kwargs["bounds"] = (lower_bounds, upper_bounds)

popt, pcov = curve_fit(ll4_formula, x_data, y_data, **fit_kwargs)
```

**Proof it isn't just theory — same optimum, from three different algorithms.**
Fitting `ryegrass` and `lepidium` (real datasets, `drc`'s own worked examples) through
R's `drm()` and through unbounded `fit_ll4` side by side, at full precision:

| dataset | param | R (BFGS) | Python (LM) | diff |
|---|---|---|---|---|
| ryegrass | b | 2.9822190713 | 2.9822229910 | 4e-6 |
| | c | 0.4814131884 | 0.4814078385 | -5e-6 |
| | d | 7.7929582937 | 7.7929624952 | 4e-6 |
| | e | 3.0579549665 | 3.0579571875 | 2e-6 |
| lepidium | b | 1.928571953 | 1.9275703918 | -0.0010 |
| | c | 3.360528391 | 3.3590041865 | -0.0015 |
| | d | 10.785121241 | 10.7854076158 | 0.0003 |
| | e | 3.657665669 | 3.6585361316 | 0.0009 |

Both are effectively the same answer. This is expected, not a coincidence: for a smooth,
well-conditioned objective with a good starting point, generic and specialized
optimizers converge to the same minimum — they just take different paths to get there.
Algorithm choice actually starts to matter when the starting point is bad or the
objective surface is badly shaped, not when the problem is this clean.

**Proof of efficiency: function evaluations to converge, same starting guess fed to
all three.**

| | ryegrass | lepidium |
|---|---|---|
| R (BFGS) | 35 function / 9 gradient evals | 30 function / 10 gradient evals |
| scipy `lm` | 26 evals | 31 evals |
| scipy `trf` | 6 evals | 7 evals |

(Confirmed `lm` and `trf` land on the same optimum as each other too — agreement to
6-7 decimal places — so `trf`'s far lower count isn't cutting corners here, it's
genuinely converging faster once already close to the answer.) These counts aren't
perfectly apples-to-apples — BFGS's numerically-differentiated gradient evaluations are
more expensive per call than a single LM/TRF function evaluation — but the gap is large
enough to be directionally real, not noise.

**Proof outlier robustness isn't about bounds at all.** A moderate outlier (one point
at 2.5x the normal max) pulls R's BFGS fit and Python's unbounded LM fit by the *same*
amount, in the same direction, to within 5e-5 on every parameter. An extreme outlier
(one point at 60x the normal max) makes *both* blow up into physically nonsensical fits
(`d` far outside the observed range, `e` collapsing near zero) — different exact
numbers, same qualitative failure. The reason: ordinary least-squares itself isn't
robust to outliers, regardless of which optimizer or which initial-guess strategy is
used — every point, including the outlier, still pulls on the final answer. Percentile-
based `c`/`d` only protects the *initial guess*; it was never a defense for the final
fit, and R has this exact same blind spot since it isn't bounded by default either.

**Literature context, not just this repo's own tests:** the tradeoff between LM and
trust-region methods generally is a documented one, not something specific to this
dataset. One direct comparison found LM converging in far fewer iterations than a
trust-region method on a given problem, but landing on a measurably *less accurate*
answer — LM's speed comes from a cheaper curvature approximation that can cut corners a
trust-region method doesn't. A separate applied comparison (battery-model parameter
fitting) found trust-region and BFGS *outperforming* LM on accuracy for that problem.
In other words, "LM is the fast one" doesn't generalize to "LM is the best one" — it
depends on the specific problem's conditioning, same as everything else in this
document.

**Verdict:** there is no universal winner between these three, and R's choice of BFGS
isn't evidence that it's the statistically superior pick for this specific job — it's
the algorithm general enough to serve `drm()`'s much broader scope (every response
distribution it supports), not one chosen for being optimal at continuous 4PL fitting
specifically. For a pure continuous least-squares problem like this one, scipy's `lm`
and `trf` are the specialized tools actually built for the job, and they prove it by
landing on R's own answer while doing meaningfully less work to get there. Between the
two: `trf` is the more robust general default (handles bounds, gives a usable
covariance matrix under rank-deficiency, and was the most efficient of all three in
testing here), while `lm` remains a reasonable choice when bounds aren't needed and the
fit will be inspected by eye regardless.
