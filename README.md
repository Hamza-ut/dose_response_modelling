**STEP 1: Correction**

There are two ways to clean up your raw OD numbers:

- **Blank Correction (Control Wells)** subtracts **across the plate** at one single timepoint.
- **Baseline Correction (Time-Zero)** subtracts **down one well** over time.

---

**Blank Correction (Control Wells)**

- **How it works:** Find the well with no bacteria (`uM = -1`) at a specific hour on a specific plate. Subtract its OD from every bacterial well measured on that plate at that exact same hour.
- **The Grouping:** Group by `Date` + `Plt/Replicate` + `Time_h` (Adding `Species` is harmless, but `Plt` already isolates the plate snapshot).
- **The Math:**
  $$OD_{\text{corrected}} = RawOD(\text{Sample Well at Hour } t) - RawOD(\text{Control Well at Hour } t)$$
- **When to use it:** When your plate has dedicated negative control wells (media only).

---

**Baseline Correction (Time-Zero)**

- **How it works:** Find a specific well's starting OD at Hour 0. Subtract that starting number from all future timepoints of that same physical well.
- **The Grouping:** Group by `Date` + `Plt/Replicate` + `Species` + `Well`.
- **The Math:**
  $$OD_{\text{corrected}} = RawOD(\text{Well at Hour } t) - RawOD(\text{Well at Hour } 0)$$
- **When to use it:** When you don't have negative control wells.

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
