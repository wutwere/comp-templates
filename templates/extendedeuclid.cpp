int64_t euclid(int64_t a, int64_t b, int64_t &x, int64_t &y) {
	if (!b) return x = 1, y = 0, a;
	int64_t d = euclid(b, a % b, y, x);
	return y -= a/b * x, d;
}
