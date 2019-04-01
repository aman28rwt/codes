
public class Treap {
	class node {
		int prior, value, cnt;
		boolean rev;
		node l, r;

		node() {
			prior = value = cnt = 0;
			rev = false;
			l = r = null;
		}
		
		node(node n) {
			this.prior = n.prior;
			this.value = n.value;
			this.cnt = n.cnt;
			this.rev = n.rev;
			this.l = n.l;
			this.r = n.r;
		}
	}

	int getCount(node n) {
		return n == null ? 0 : n.cnt;
	}

	void updateCount(node n) {
		if (n != null) {
			n.cnt = getCount(n.l) + getCount(n.r) + 1;
		}
	}

	void push(node n) {
		if (n != null && n.rev) {
			n.rev = false;

			node temp = n.l;
			n.l = n.r;
			n.r = temp;

			if (n.l != null) {
				n.l.rev = !n.l.rev;
			}
			if (n.r != null) {
				n.r.rev = !n.r.rev;
			}
		}
	}

	void merge(node n, node l, node r) {
		push(l);
		push(r);
		if ((l == null && r != null) || (l != null && r == null)) {
			n = l != null ? l : r;
		} else if (l.prior > r.prior) {
			merge(l.r, l.r, r);
			n = new node(l);
		} else {
			merge(r.l, l, r.l);
			n = new node(r);
		}
		updateCount(n);
	}

	void split(node n, node l, node r, int key, int add) {
		if (n == null) {
			l = r = null;
		}
		push(n);
		int curr_key = add + getCount(n.l);
		if (key < curr_key) {
			split(n.l, l, n.l, key, add);
			r = new node(n);
		} else {
			split(n.r, n.r, r, key, add + 1 + getCount(n.l));
			l = new node(n);
		}
		updateCount(n);
	}

	void reverse(node n, int l, int r) {
		node n1 = new node(), n2 = new node(), n3 = new node();
		split(n, n1, n2, l, 0);
		split(n2, n2, n3, r - l + 1, 0);
		n2.rev = !n2.rev;
		merge(n, n1, n2);
		merge(n, n, n3);
	}

	String print(node n) {
		if (n == null) {
			return "";
		}
		push(n);
		return print(n.l) + " " + n.value + " " + print(n.r);
	}
}