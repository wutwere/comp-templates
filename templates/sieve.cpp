const int MAXN = 50000000;
bitset<MAXN> isprime;

isprime.set();
isprime[0] = isprime[1] = 0;
for (int i = 4; i < MAXN; i += 2) isprime[i] = 0;
for (int i = 3; i * i < MAXN; i += 2)
  if (isprime[i])
    for (int j = i * i; j < MAXN; j += i * 2)
      isprime[j] = 0;
