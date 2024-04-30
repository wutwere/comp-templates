const int n = 1e6;
int64_t phi[n + 1];

iota(phi, phi + n + 1, -1);
phi[1] = 1;
for (int i = 2; i <= n; i++)
  if (phi[i] == i - 1)
    for (int j = i + i; j <= n; j += i)
      phi[j] -= phi[j] / i;
