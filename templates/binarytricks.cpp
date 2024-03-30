// Print x in binary
void print(int x) {
  cout << bitset<32>(x) << '\n';
}

int x = 1234;
// Get least set bit in x
x & -x;

// Iterate over all subsets of set bits in x, except x itself
for (int s = x; s; ) {
  --s &= x;
  // s is the subset
}

// Get the next number after x that has the same number of set bits
int c = x & -x;
int r = x + c;
r |= ((r ^ x) >> 2) / c;

// For all i up to (1 << M), compute the sum of all subsets of set bits of i
int M = 5;
vector<int> d(1 << M);
iota(d.begin(), d.end(), 0);
for (int b = 0; b < M; b++)
  for (int i = 0; i < (1 << M); i++)
    if (i >> b & 1)
      d[i] += d[i ^ (1 << b)];
