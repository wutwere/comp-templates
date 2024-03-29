// n choose r mod MOD in O(1) after O(n) precomputation
int64_t fac[MAXN + 1], finv[MAXN + 1], inv[MAXN + 1];
int64_t ncr(int n, int r) {
  if (fac[0] == 0) {
    fac[0] = finv[0] = inv[1] = 1;
    for (int i = 2; i <= MAXN; i++) {
      inv[i] = inv[MOD % i] * (-MOD / i) % MOD;
      if (inv[i] < 0) inv[i] += MOD;
    }
    for (int i = 1; i <= MAXN; i++) {
      fac[i] = fac[i - 1] * i % MOD;
      finv[i] = finv[i - 1] * inv[i] % MOD;
    }
  }
  return (n < 0 || r < 0 || r > n ? 0 : fac[n] * finv[r] % MOD * finv[n - r] % MOD);
}
