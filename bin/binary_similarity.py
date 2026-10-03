
import numpy as np
from typing import Iterable, Dict, Any, Tuple

def _to_bool_array(x: Iterable) -> np.ndarray:
    """Convert an iterable to a 1D boolean numpy array of dtype=bool."""
    arr = np.asarray(x).astype(bool).ravel()
    return arr

def _confusion_counts(a: np.ndarray, b: np.ndarray) -> Tuple[int, int, int, int]:
    """Return TP, FP, FN, TN for two boolean arrays."""
    tp = int(np.sum(a & b))
    tn = int(np.sum(~a & ~b))
    fp = int(np.sum(a & ~b))
    fn = int(np.sum(~a & b))
    return tp, fp, fn, tn

def _entropy_dist(p: np.ndarray) -> float:
    p = np.asarray(p, dtype=float)
    p = p[p > 0.0]
    if p.size == 0:
        return 0.0
    return float(-np.sum(p * np.log2(p)))

def _kl_bernoulli(p: float, q: float) -> float:
    def term(a, b):
        if a == 0.0:
            return 0.0
        if b == 0.0:
            return np.inf
        return a * (np.log2(a) - np.log2(b))
    return float(term(p, q) + term(1.0 - p, 1.0 - q))

def _jsd_bernoulli(p: float, q: float) -> float:
    m = 0.5 * (p + q)
    return 0.5 * _kl_bernoulli(p, m) + 0.5 * _kl_bernoulli(q, m)

def compute_binary_similarity(a: Iterable, b: Iterable) -> Dict[str, Any]:
    A = _to_bool_array(a)
    B = _to_bool_array(b)
    if A.shape != B.shape:
        raise ValueError(f"Vectors must have the same shape, got {A.shape} and {B.shape}")
    n = A.size
    tp, fp, fn, tn = _confusion_counts(A, B)
    def safe_div(num, den):
        return float(num) if den == 0 else float(num) / float(den)
    dice = safe_div(2 * tp, 2 * tp + fp + fn)
    jaccard = safe_div(tp, tp + fp + fn)
    overlap = safe_div(tp, min(tp + fp, tp + fn))
    hamming_similarity = safe_div(tp + tn, n)
    denom = np.sqrt((tp + fp) * (tp + fn))
    cosine = safe_div(tp, denom)
    denom_phi = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    pearson_phi = safe_div(tp * tn - fp * fn, denom_phi)
    p11 = tp / n
    p10 = fp / n
    p01 = fn / n
    p00 = tn / n
    px1 = p11 + p10
    py1 = p11 + p01
    Hx = _entropy_dist(np.array([px1, 1.0 - px1]))
    Hy = _entropy_dist(np.array([py1, 1.0 - py1]))
    MI = 0.0
    joint = np.array([[p00, p01],[p10, p11]], dtype=float)
    px = np.array([1.0 - px1, px1], dtype=float)
    py = np.array([1.0 - py1, py1], dtype=float)
    for ix in (0,1):
        for iy in (0,1):
            pxy = joint[ix, iy]
            if pxy > 0.0:
                MI += pxy * (np.log2(pxy) - np.log2(px[ix]) - np.log2(py[iy]))
    NMI_sqrt = (MI / np.sqrt(Hx * Hy)) if (Hx > 0 and Hy > 0) else 0.0
    NMI_avg = (2.0 * MI / (Hx + Hy)) if (Hx + Hy) > 0 else 0.0
    VI = Hx + Hy - 2.0 * MI
    KL_XY = _kl_bernoulli(px1, py1)
    KL_YX = _kl_bernoulli(py1, px1)
    JSD = _jsd_bernoulli(px1, py1)
    return {
        "n": n, "tp": tp, "fp": fp, "fn": fn, "tn": tn,
        "dice": dice, "jaccard": jaccard, "overlap": overlap,
        "hamming_similarity": hamming_similarity, "cosine": cosine, "pearson_phi": pearson_phi,
        "Hx_bits": Hx, "Hy_bits": Hy, "MI_bits": MI, "NMI_sqrt": NMI_sqrt, "NMI_avg": NMI_avg, "VI_bits": VI,
        "KL_XY_bits": KL_XY, "KL_YX_bits": KL_YX, "JSD_bits": JSD,
    }

def compute_info_measures(a: Iterable, b: Iterable) -> Dict[str, Any]:
    res = compute_binary_similarity(a, b)
    keys = ["Hx_bits","Hy_bits","MI_bits","NMI_sqrt","NMI_avg","VI_bits","KL_XY_bits","KL_YX_bits","JSD_bits"]
    return {k: res[k] for k in keys}
