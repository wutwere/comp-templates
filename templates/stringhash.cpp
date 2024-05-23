string s, t; cin >> s >> t;
const int64_t M = 1e9 + 9;
mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
int64_t B = uniform_int_distribution<int64_t>(0, M - 1)(rng);
vector<int64_t> pow{1};
for (int i = 0; i < s.size(); i++) {
  pow.push_back((pow.back() * B) % M);
}
auto hashes = [&](string &s) -> vector<int64_t> {
  vector<int64_t> h(s.size() + 1);
  for (int i = 0; i < s.size(); i++) {
    h[i + 1] = (h[i] * B + s[i]) % M;
  }
  return h;
};
auto hash_range = [&](vector<int64_t> &hashes, int l, int r) -> int64_t {
  // half open range [l, r)
  int64_t h = hashes[r] - hashes[l] * pow[r - l] % M;
  if (h < 0) h += M;
  return h;
};
auto h1 = hashes(s), h2 = hashes(t);
for (int i = 0; i + 2 <= s.size(); i++) {
  if (hash_range(h1, i, i + t.size()) == h2.back()) {
    cout << i << '\n';
  }
}
