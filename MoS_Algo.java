public class MoS_Algo {
	
	public static int blockSize;
	
	public static class query implements Comparable<query> {
		int l, r, i;

		@Override
		public int compareTo(query o) {
			int mul = (this.l / blockSize) % 2 == 0 ? 1 : -1;
			return this.l / blockSize == o.l / blockSize ? (this.r - o.r) * mul : this.l / blockSize - o.l / blockSize;
		}
	}
	
	
}