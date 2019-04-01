import java.util.Arrays;

public class LongHashCounter {
	public long[] keys;
	public int[] allocated;
	private int scale = 1 << 2;
	private int rscale = 1 << 1;
	private int mask = scale - 1;
	public int size = 0;

	public LongHashCounter() {
		allocated = new int[scale];
		Arrays.fill(allocated, NG);
		keys = new long[scale];
	}

	// if value is NG, entry is removed. (e.g. 0)
	private static final int NG = 0;

	public boolean containsKey(long x) {
		int pos = h(x) & mask;
		while (allocated[pos] != NG) {
			if (x == keys[pos])
				return true;
			pos = pos + 1 & mask;
		}
		return false;
	}

	public int get(long x) {
		int pos = h(x) & mask;
		while (allocated[pos] != NG) {
			if (x == keys[pos])
				return allocated[pos];
			pos = pos + 1 & mask;
		}
		return NG;
	}

	public int put(long x, int v) {
		int pos = h(x) & mask;
		while (allocated[pos] != NG) {
			if (x == keys[pos]) {
				int oldval = allocated[pos];
				allocated[pos] = v;
				return oldval;
			}
			pos = pos + 1 & mask;
		}
		if (size == rscale) {
			resizeAndPut(x, v);
		} else {
			keys[pos] = x;
			allocated[pos] = v;
		}
		size++;
		return NG;
	}

	public int inc(long x, int v) {
		int pos = h(x) & mask;
		while (allocated[pos] != NG) {
			if (x == keys[pos]) {
				allocated[pos] += v;
				return allocated[pos];
			}
			pos = pos + 1 & mask;
		}
		if (size == rscale) {
			resizeAndPut(x, v);
		} else {
			keys[pos] = x;
			allocated[pos] = v;
		}
		size++;
		return v;
	}

	public boolean remove(long x) {
		int pos = h(x) & mask;
		while (allocated[pos] != NG) {
			if (x == keys[pos]) {
				size--;
				// take last and fill rmpos
				int last = pos;
				pos = pos + 1 & mask;
				while (allocated[pos] != NG) {
					int lh = h(keys[pos]) & mask;
					// lh <= last < pos
					if (lh <= last && last < pos || pos < lh && lh <= last || last < pos && pos < lh) {
						keys[last] = keys[pos];
						allocated[last] = allocated[pos];
						last = pos;
					}
					pos = pos + 1 & mask;
				}
				keys[last] = 0;
				allocated[last] = NG;

				return true;
			}
			pos = pos + 1 & mask;
		}
		return false;
	}

	private void resizeAndPut(long x, int v) {
		int nscale = scale << 1;
		int nrscale = rscale << 1;
		int nmask = nscale - 1;
		int[] nallocated = new int[nscale];
		Arrays.fill(nallocated, NG);
		long[] nkeys = new long[nscale];
		for (int i = next(0); i < scale; i = next(i + 1)) {
			long y = keys[i];
			int pos = h(y) & nmask;
			while (nallocated[pos] != NG)
				pos = pos + 1 & nmask;
			nkeys[pos] = y;
			nallocated[pos] = allocated[i];
		}
		{
			int pos = h(x) & nmask;
			while (nallocated[pos] != NG)
				pos = pos + 1 & nmask;
			nkeys[pos] = x;
			nallocated[pos] = v;
		}
		allocated = nallocated;
		keys = nkeys;
		scale = nscale;
		rscale = nrscale;
		mask = nmask;
	}

	public int next(int itr) {
		while (itr < scale && allocated[itr] == NG)
			itr++;
		return itr;
	}

	private int h(long x) {
		x ^= x >>> 33;
		x *= 0xff51afd7ed558ccdL;
		x ^= x >>> 33;
		x *= 0xc4ceb9fe1a85ec53L;
		x ^= x >>> 33;
		return (int) x;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (int i = next(0); i < scale; i = next(i + 1)) {
			sb.append(",");
			sb.append(keys[i] + ":" + allocated[i]);
		}
		return sb.length() == 0 ? "" : sb.substring(1);
	}
}