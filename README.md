**STEP 1: Correction**

There are two ways to clean up your raw OD numbers:

- **Blank Correction (Media/Control Wells)** subtracts **across the plate** at one single timepoint.
- **Baseline Correction (Minimum Value)** subtracts **down one well** over time, using that well's own lowest reading.

---

**Blank Correction (Media/Control Wells)**

- **How it works:** Find the well with no bacteria (`uM = -1`) on every specific plate/replicate for each hour. Subtract its OD from every bacterial well measured on that plate at that exact same hour.
- **The Grouping:** Group by `Plt/Replicate` + `Species` + `Time_h`.
- **The Math:**
  $$OD_{\text{corrected}} = RawOD(\text{Sample Well at Hour } t) - RawOD(\text{Control Well at Hour } t)$$
- **When to use it:** When your data has dedicated media only wells.

---

**Baseline Correction (Minimum Value)**

- **How it works:** For each unique Replicate & Well, find LOWEST OD value across its entire time series (not necessarily at earliest time/ hour 0 of that replicate & well) . Subtract that minimum from every timepoint of that same physical well. This is more robust than using the Hour-0 reading alone, since any single reading can be a noisy outlier (a bubble, a pipetting blip) -- using the well's own lowest point means one bad reading can't set a bad baseline for the whole series.
- **The Grouping:** Group by whatever uniquely identifies one physical well -- e.g. `Date` + `Plt/Replicate` + `Species` + `Well` (or `Plt_well`). Verify this actually isolates one well for your dataset: `df.groupby(group_fields + [dose_field])[well_field].nunique().max()` should be `1` if `dose_field` alone is enough; if it's `>1`, real well identity is required in the grouping.
- **The Math:**
  $$OD_{\text{corrected}} = RawOD(\text{Well at Hour } t) - \min_t\big(RawOD(\text{Well})\big)$$
- **When to use it:** When you don't have negative control wells, and you don't expect background to drift over time. Matches the default `min` method in `growthcurver`, a widely-used R package for bacterial growth curve analysis.

---

---

**STEP 2: Biological Aggregation/Normalization (`OD_normalized`)**

After Step 1, you get the OD_corrected values (either way1 or way2 depends on ur dataset) which are cleaned of background noise. Step 2 converts those OD_corrected growth into a normalized scale (**0.0 to 1.0**) for each group/replicate.

To measure drug effectiveness, you must compare the drug-treated well against a baseline of **100% healthy, uninhibited growth** (`uM = 0`):

$$\text{OD}_{\text{normalized}} = \frac{\text{OD}_{\text{corrected}}}{\text{OD}_{\text{corrected}}(\text{uM} = 0)}$$

We use uM=0 of the corresponding plate to calculate `OD_normalized`. Grouping can be as below

1. **`Date`** (Experiment run(if any))
2. **`Plt / Replicate`** (Physical plate ID)
3. **`Time_h`** (Exact hour of measurement)

> **The Rule:** You ONLY compare a treated well against the `uM = 0` control well on its **OWN physical plate** at its **OWN hour**.

**What the Output Numbers Mean**

- **`1.0` (100%):** Equal to normal growth — the drug had zero effect.
- **`0.5` (50%):** Half growth — the drug reduced growth by 50%.
- **`0.0` (0%):** Zero growth — the drug completely stopped/killed the bacteria.

---

---

**Step 3: Modelling 4PL on your desired timepoint**

---

---

**Step 4:**
