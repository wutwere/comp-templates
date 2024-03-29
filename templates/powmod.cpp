// x to the y-th power mod MOD in O(log y)
int powmod(int x, int y) {
  int res = 1;
  for (; y; y >>= 1, x = 1ll * x * x % MOD)
    if (y & 1) res = 1ll * res * x % MOD;
  return res % MOD;
}
