"""BH and Bonferroni adjustment for a supplied p-value family."""
import numpy as np

def adjust_global(pvalues):
    """Return BH adjusted p-values and Bonferroni adjusted p-values."""
    values = np.asarray(pvalues, dtype=np.float64)
    if values.ndim != 1 or not len(values):
        raise ValueError("Expected a nonempty vector of p-values")
    if not np.isfinite(values).all() or ((values < 0) | (values > 1)).any():
        raise ValueError("P-values must be finite and within [0, 1]")
    m = len(values)
    order = np.argsort(values, kind="stable")
    ranked = values[order]
    adjusted_sorted = np.minimum.accumulate(
        (ranked * m / np.arange(1, m + 1))[::-1])[::-1]
    adjusted_sorted = np.minimum(adjusted_sorted, 1.0)
    bh = np.empty(m, dtype=np.float64)
    bh[order] = adjusted_sorted
    bonferroni = np.minimum(values * m, 1.0)
    if (np.diff(adjusted_sorted) < -1e-14).any():
        raise AssertionError("BH adjusted values are not monotone")
    if (bh < values - 1e-12).any() or (bonferroni < bh - 1e-12).any():
        raise AssertionError("Invalid ordering of adjusted p-values")
    return bh, bonferroni
