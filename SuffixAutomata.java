import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;

public class SuffixAutomata {

	public static class Suffix {
		public static class node {
			int len, link;
			HashMap<Character, Integer> next = new HashMap<>();
		}

		node[] arr;
		int size, last;

		Suffix(int len) {
			this.arr = new node[2 * len + 1];
			for (int i = 0; i < arr.length; i++) {
				arr[i] = new node();
			}
			this.arr[0].len = 0;
			this.arr[0].link = -1;
			this.size = 1;
			this.last = 0;
		}

		public void extend_SA(Character ch) {
			int curr = size++;
			arr[curr].len = arr[last].len + 1;
			int p;
			for (p = last; p != -1 && !arr[p].next.containsKey(ch); p = arr[p].link) {
				arr[p].next.put(ch, curr);
			}

			if (p == -1) {
				arr[curr].link = 0;
			} else {
				int q = arr[p].next.get(ch);
				if (arr[p].len + 1 == arr[q].len) {
					arr[curr].link = q;
				} else {
					int clone = size++;
					arr[clone].len = arr[p].len + 1;
					arr[clone].next = arr[q].next;
					arr[clone].link = arr[q].link;

					for (; p != -1 && arr[p].next.get(ch) == q; p = arr[p].link) {
						arr[p].next.put(ch, clone);
					}
					arr[q].link = arr[curr].link = clone;
				}
			}

			last = curr;
		}
	}

	public static ArrayList<Character> list = new ArrayList<>();
	public static int path = 0;
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		Scanner scn = new Scanner(System.in);

		System.out.println("GO");

		String str = scn.next();
		Suffix sf = new Suffix(str.length());

		for (int i = 0; i < str.length(); i++) {
			sf.extend_SA(str.charAt(i));
		}

		kLex(0, scn.nextInt() - 1, sf);
		System.out.println(list);
	}

	public static void kLex(int node, int k, Suffix sf) {
		for (int i = 0; i < 26; i++) {
			if (sf.arr[node].next.containsKey((char)(i + 'a'))) {
				path++;
				if (path == k) {
					list.add(0, (char) ('a' + i));
					return;
				}
				kLex(sf.arr[node].next.get((char) (i + 'a')), k, sf);
				if (path == k) {
					list.add(0, (char) ('a' + i));
					return;
				}
			}
		}
	}
}