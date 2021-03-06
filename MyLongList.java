import java.util.Arrays;
import java.util.Random;

public class MyLongList {
	private long[] arr;
	private int capacity = 1 << 4;
	private int initialCapacity = 1 << 4;
	public int length = 0;
	private int lambda = 2;

	public MyLongList() {
		arr = new long[capacity];
	}
	
	public MyLongList(int capacity) {
		this.initialCapacity = this.capacity = capacity;
		arr = new long[capacity];
	}
	
	public MyLongList(long[] arr) {
		this.length = arr.length;
		this.capacity = lambda * arr.length;
		this.arr = arr;
	}
	
	public void clear() {
		capacity = initialCapacity;
		arr = new long[capacity];
	}
	
	public boolean isEmpty() {
		return length == 0;
	}
	
	public int size() {
		return this.length;
	}
	
	public int capacity() {
		return this.capacity;
	}
	
	public void add(long val) {
		if(length == capacity) {
			grow();
		}
		arr[length++] = val;
	}
	
	public void removeLast() {
		arr[--length] = 0;
		if(capacity / lambda >= initialCapacity && length < capacity / (1.5 * lambda)) {
			trim();
		}
	}
	
	private void trim() {
		capacity /= lambda;
		arr = Arrays.copyOf(arr, capacity);
	}
	
	private void grow() {
		capacity *= lambda;
		arr = Arrays.copyOf(arr, capacity);
	}
	
	public long get(int ind) {
		if(ind < 0 || ind >= length) {
			throw new ArrayIndexOutOfBoundsException();
		}
		return arr[ind];
	}
	
	public void set(int ind, long val) {
		if(ind < 0 || ind >= length) {
			throw new ArrayIndexOutOfBoundsException();
		}
		arr[ind] = val;
	}
	
	public void sort() {
		Random r = new Random();
		for (int i = 1, j; i < length; i++) {
			j = r.nextInt(i);
			arr[i] = arr[i] ^ arr[j];
			arr[j] = arr[i] ^ arr[j];
			arr[i] = arr[i] ^ arr[j];
		}
		Arrays.sort(arr, 0, length);
	}
}
