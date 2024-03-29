// vector<int> m = manacher(s);
// m[2 * i] for even length palindrome centered between i - 1 and i
// m[2 * i + 1] for odd length palindrome centered at i
vector<int> manacher(string &s) {
  string odd = "#";
  for (char c : s) odd += c, odd += '#';
  odd = "$" + odd + "^";
  vector<int> m(odd.size());
  for (int i = 1, l = 1, r = 1; i < odd.size(); i++) {
    if (i < r) m[i] = max(0, min(r - i, m[l + r - i]));
    if (i + m[i] >= r) while (odd[i - m[i] - 1] == odd[i + m[i] + 1])
      m[i]++, l = i - m[i], r = i + m[i];
  }
  return vector<int>(m.begin() + 1, m.end() - 2);
}
