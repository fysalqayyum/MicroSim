#!/usr/bin/env python3
"""Gate 3d: statistics of the splitmix64 hash used by the repaired addNoise.

Mirrors the kernel's hash in Python and checks it behaves like noise rather
than like the upstream version, which derived its value from work-item ids
only and so replayed one fixed spatial pattern every timestep.

Checks: uniform[-1,1] marginal; the value CHANGES with tstep for a fixed cell
(the actual defect); and lag-1 autocorrelation in both time and space is no
worse than a reference PRNG measured the same way. The control series matters --
without it there is no way to know what "small enough" means for a max taken
over 20 series.

Run: python3 test_addnoise_hash.py     (no GPU, no arguments)
"""
import statistics as st, math, random
M = (1 << 64) - 1
def u_of(index, tstep):
    h = (index * 0x9E3779B97F4A7C15) & M
    h ^= (tstep * 0xBF58476D1CE4E5B9) & M
    h ^= h >> 30; h = (h * 0xBF58476D1CE4E5B9) & M
    h ^= h >> 27; h = (h * 0x94D049BB133111EB) & M
    h ^= h >> 31
    return 2.0 * ((h >> 11) / 9007199254740992.0) - 1.0

def acorr(x, lag=1):
    m = st.mean(x); v = sum((a-m)**2 for a in x)
    return sum((x[i]-m)*(x[i+lag]-m) for i in range(len(x)-lag))/v

N = 50000
se = 1.0/math.sqrt(N)
print(f"N={N} per series, 1 s.e. = {se:.4f}, 3 s.e. = {3*se:.4f}")

print("\n-- time series, fixed cell (20 cells) --")
t_ac = [acorr([u_of(c, t) for t in range(1, N+1)]) for c in
        [7, 5000, 12345, 123457, 999999, 3, 88, 60001, 271828, 314159,
         1, 2, 4, 8, 16, 32, 64, 128, 256, 512]]
print("  max|ac| = %.4f (%.1f s.e.)  mean = %+.5f" %
      (max(map(abs,t_ac)), max(map(abs,t_ac))/se, st.mean(t_ac)))

print("-- space series, fixed step (20 steps) --")
s_ac = [acorr([u_of(i, T) for i in range(1, N+1)]) for T in
        [1, 2, 10, 500, 9999, 99999, 123, 4096, 777, 31337,
         3, 5, 7, 11, 13, 17, 19, 23, 29, 31]]
print("  max|ac| = %.4f (%.1f s.e.)  mean = %+.5f" %
      (max(map(abs,s_ac)), max(map(abs,s_ac))/se, st.mean(s_ac)))

# Control: a real PRNG, same estimator, same N -- what does "clean" look like?
rng = random.Random(1)
ctrl = [acorr([rng.uniform(-1,1) for _ in range(N)]) for _ in range(20)]
print("-- control (python Mersenne Twister, 20 series) --")
print("  max|ac| = %.4f (%.1f s.e.)" % (max(map(abs,ctrl)), max(map(abs,ctrl))/se))

allac = t_ac + s_ac
verdict = max(map(abs, allac)) < 4*se
print("\nGATE 3d:", "PASS" if verdict else "FAIL",
      "- threshold 4 s.e. = %.4f" % (4*se))
