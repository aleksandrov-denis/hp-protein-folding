// From http://www.cburch.com/ccsc-nifty/12/bur/handout.html

import java.io.PrintStream;
import java.util.Collections;
import java.util.Arrays;
import java.lang.Math;
import java.util.concurrent.ThreadLocalRandom;
import java.util.ArrayList;

public class Protein implements Cloneable {
	public String chain; // string of H or P's representing amino acids
	public int[] row;    // row[i] represents y-coordinate of residue i
	public int[] col;    // col[i] represents x-coordinate of residue i
	public ArrayList<Integer> corners, crankshafts, flats;
	static int best_score;

	/** Constructs a protein that is simply lined. */
	public Protein(String chain) {
		for (int i = 0; i < chain.length(); i++) {
			char c = chain.charAt(i);
			if (c != 'H' && c != 'P') {
				throw new IllegalArgumentException("must all be P or H");
			}
		}
		this.chain = chain;
		this.row = new int[chain.length()];
		this.col = new int[chain.length()];
		for (int i = 0; i < col.length; i++) {
			row[i] = 0;
			col[i] = i;
		}
		corners = new ArrayList<Integer>();
		crankshafts = new ArrayList<Integer>();
		flats = new ArrayList<Integer>();
	}

	/** Returns the number of amino acids in this protein chain. */
	public int getLength() {
		return chain.length();
	}

	/** Returns true if residue at index is hydrophobic. */
	public boolean isHydrophobic(int index) {
		return chain.charAt(index) == 'H';
	}

	/** Returns y-coordinate at which residue at index has been placed. */
	public int getRow(int index) {
		return row[index];
	}

	/** Returns x-coordinate at which residue at index has been placed. */
	public int getColumn(int index) {
		return col[index];
	}

	/** Returns the height of the protein **/
	private int[] getHeight(int start) {
		int[] out = new int[2];
		int max = Integer.MIN_VALUE;
		int min = Integer.MAX_VALUE;
		for (int i = start; i < this.chain.length(); i++) {
			if (this.row[i] < min) {
				min = this.row[i];
			}
			if (this.row[i] > max) {
				max = this.row[i];
			}
		}
		out[0] = max;
		out[1] = min;
		return out;
	}

	/** Returns the width of the protein **/
	private int[] getWidth(int start) {
		int[] out = new int[2];
		int max = Integer.MIN_VALUE;
		int min = Integer.MAX_VALUE;
		for (int i = start; i < this.chain.length(); i++) {
			if (this.col[i] < min) {
				min = this.col[i];
			}
			if (this.col[i] > max) {
				max = this.col[i];
			}
		}
		out[0] = max;
		out[1] = min;
		return out;
	}

	/** Pivot move 90 degrees clockwise **/
	public void pivot_move() {
		int k = ThreadLocalRandom.current().nextInt(0, this.chain.length());
		int[] max_min = getHeight(k);
		int[] row2 = new int[chain.length() - k];
		int[] col2 = new int[chain.length() - k];
		for (int i = k; i < chain.length(); i++) {
			row2[i - k] = this.col[i];
			col2[i - k] = (max_min[0] - max_min[1]) - (this.row[i] - max_min[1]);
		}
		for (int i = k; i < chain.length(); i++) {
			this.row[i] = row2[i - k] + (this.row[k] - row2[0]);
			this.col[i] = col2[i - k] + max_min[1] + (this.col[k] - (col2[0] + max_min[1]));
		}

	}

	/** Find all corners in current structure **/
	public void establish_corners() {
		this.corners = new ArrayList<Integer>();
		for (int k = 1; k < this.chain.length() - 1; k++) {
			if (this.row[k - 1] != this.row[k + 1] && this.col[k - 1] != this.col[k + 1]) {
				corners.add(k);
			}
		}
	}

	/** Get corner conditions **/
	public boolean[] getConditions(int r, int c, int corner) {
		boolean[] conditions = new boolean[4];
		int r_prev = this.row[corner - 1];
		int c_prev = this.col[corner - 1];
		int r_next = this.row[corner + 1];
		int c_next = this.col[corner + 1];
		conditions[0] = ((r_prev == r - 1) && (c_prev == c) && (r_next == r) && (c_next == c + 1)) || ((r_prev == r) && (c_prev == c + 1) && (r_next == r - 1) && (c_next == c));
		conditions[1] = ((r_prev == r) && (c_prev == c - 1) && (r_next == r - 1) && (c_next == c)) || ((r_prev == r - 1) && (c_prev == c) && (r_next == r) && (c_next == c - 1));
		conditions[2] = ((r_prev == r + 1) && (c_prev == c) && (r_next == r) && (c_next == c + 1)) || ((r_prev == r) && (c_prev == c + 1) && (r_next == r + 1) && (c_next == c));
		conditions[3] = ((r_prev == r) && (c_prev == c - 1) && (r_next == r + 1) && (c_next == c)) || ((r_prev == r + 1) && (c_prev == c) && (r_next == r) && (c_next == c - 1));
		return conditions;
	}

	/** Corner perturbation **/
	public void perterb_corner() {
		if (this.corners.size() == 0) {
			return;
		}
		int index = ThreadLocalRandom.current().nextInt(0, this.corners.size());
		int corner = this.corners.get(index);
		int r = this.row[corner];
		int c = this.col[corner];
		boolean[] conditions = getConditions(r, c, corner);
		if (conditions[0]) {
			// perturbation #1: bottom-left corner
			this.row[corner] = r - 1;
			this.col[corner] = c + 1;
		} else if (conditions[1]) {
			// perturbation #2: bottom-right corner
			this.row[corner] = r - 1;
			this.col[corner] = c - 1;
		} else if (conditions[2]) {
			// perturbation #3: top-left corner
			this.row[corner] = r + 1;
			this.col[corner] = c + 1;
		} else if (conditions[3]) {
			// perturbation #4: top-right corner
			this.row[corner] = r + 1;
			this.col[corner] = c - 1;
		}
	}

	/** Determine where the flat parts are **/
	public void establish_flats() {
		this.flats = new ArrayList<Integer>();
		for (int k = 2; k < this.chain.length() - 3; k++) {
			// if all six residues are in a straight line
			if ((this.row[k - 2] == this.row[k - 1] && this.row[k - 1] == this.row[k] && this.row[k] == this.row[k + 1] && this.row[k + 1] == this.row[k + 2] && this.row[k + 2] == this.row[k + 3])
					|| (this.col[k - 2] == this.col[k - 1] && this.col[k - 1] == this.col[k] && this.col[k] == this.col[k + 1] && this.col[k + 1] == this.col[k + 2] && this.col[k + 2] == this.col[k + 3])) {
				flats.add(k);
			}
		}
	}

	/** Create a new crankshaft **/
	public void make_crankshaft() {
		int direction = 0;
		if (this.flats.size() == 0) {
			return;
		}
		int index = ThreadLocalRandom.current().nextInt(0, this.flats.size());
		int k = this.flats.get(index);
		Protein p = (Protein) clone();
		do {
			direction = ThreadLocalRandom.current().nextInt(-1, 2);
		} while (direction == 0);
		// check if vertical or horizontal
		if (p.row[k] == p.row[k - 1]) {
			// horizontal
			p.row[k] += direction;
			p.row[k + 1] += direction;
			// need to reconnect <k and >k+1
			if (p.col[k - 1] < p.col[k]) {
				direction = -1;
			} else {
				direction = 1;
			}
			for (int i = k; i < this.chain.length(); i++) {
				if (i > k + 1) {
					p.col[i] += direction;
				}
				p.col[i] += direction;
			}
		} else {
			// vertical
			p.col[k] += direction;
			p.col[k + 1] += direction;
			// need to reconnect <k and >k+1
			if (p.row[k - 1] < p.row[k]) {
				direction = -1;
			} else {
				direction = 1;
			}
			for (int i = k; i < this.chain.length(); i++) {
				if (i > k + 1) {
					p.row[i] += direction;
				}
				p.row[i] += direction;
			}
		}

		if (p.getScore() != -1) {
			copy_arrays(p);
		}

	}

	/** Determine where the crankshafts are **/
	public void establish_crankshafts() {
		this.crankshafts = new ArrayList<Integer>();
		for (int k = 2; k < this.chain.length() - 3; k++) {
			// if (k - 1) amd (k + 2) are corners
			if (this.row[k - 2] != this.row[k] && this.col[k - 2] != this.col[k] && this.row[k + 1] != this.row[k + 3] && this.col[k + 1] != this.col[k + 3]) {
				// sufficient to check if (k-2),(k-1),(k+2),(k+3) are colinear
				// no need to check colinearity of k,(k+1) since the structures are valid
				if (this.row[k - 2] == this.row[k - 1] && this.row[k - 2] == this.row[k + 2] && this.row[k - 2] == this.row[k + 3]) {
					crankshafts.add(k);
				} else if (this.col[k - 2] == this.col[k - 1] && this.col[k - 2] == this.col[k + 2] && this.col[k - 2] == this.col[k + 3]) {
					crankshafts.add(k);
				}
			}
		}
	}

	/** Make a crankshaft perturbation **/
	public void perterb_crank() {
		if (this.crankshafts.size() == 0) {
			return;
		}
		int index = ThreadLocalRandom.current().nextInt(0, this.crankshafts.size());
		int corner0 = this.crankshafts.get(index);
		int r0 = this.row[corner0];
		int c0 = this.col[corner0];
		boolean[] conditions0 = getConditions(r0, c0, corner0);
		int corner1 = corner0 + 1;
		int r1 = this.row[corner1];
		int c1 = this.col[corner1];
		boolean[] conditions1 = getConditions(r1, c1, corner1);
		if (conditions0[2] && conditions1[0] || conditions0[0] && conditions1[2]) {
			// perturbation a: vertical crankshaft, corners on the left
			this.row[corner0] = r0;
			this.col[corner0] = c0 + 2;
			this.row[corner1] = r1;
			this.col[corner1] = c1 + 2;
		} else if (conditions0[3] && conditions1[1] || conditions0[1] && conditions1[3]) {
			// perturbation b: vertical crankshaft, corners on the right
			this.row[corner0] = r0;
			this.col[corner0] = c0 - 2;
			this.row[corner1] = r1;
			this.col[corner1] = c1 - 2;
		} else if (conditions0[0] && conditions1[1] || conditions0[1] && conditions1[0]) {
			// perturbation c: horizontal crankshaft, corners on the bottom
			this.row[corner0] = r0 - 2;
			this.col[corner0] = c0;
			this.row[corner1] = r1 - 2;
			this.col[corner1] = c1;
		} else if (conditions0[2] && conditions1[3] || conditions0[3] && conditions1[2]) {
			// perturbation d: horizontal crankshaft, corners on the top
			this.row[corner0] = r0 + 2;
			this.col[corner0] = c0;
			this.row[corner1] = r1 + 2;
			this.col[corner1] = c1;
		}
	}

	/** Rotates the structure 90 degrees **/
	public void rotate() {
		int[] max_min = getHeight(0);
		int[] row2 = new int[chain.length()];
		int[] col2 = new int[chain.length()];
		for (int i = 0; i < chain.length(); i++) {
			row2[i] = this.col[i];
			col2[i] = (max_min[0] - max_min[1]) - (this.row[i] - max_min[1]);
		}
		for (int i = 0; i < chain.length(); i++) {
			this.row[i] = row2[i];
			this.col[i] = col2[i] + max_min[1];
		}
	}

	public void printArrays() {
		System.out.println("Row:");
		for (int i = 0; i < chain.length(); i++) {
			System.out.print(row[i] + " ");
		}
		System.out.println("\nCol:");
		for (int i = 0; i < chain.length(); i++) {
			System.out.print(col[i] + " ");
		}
		System.out.println();
	}

	/** Compares two protein structures **/
	public boolean isEqual(Protein pr) {
		int min_height = getHeight(0)[1];
		int min_width = getWidth(0)[1];
		int p_min_height = pr.getHeight(0)[1];
		int p_min_width = pr.getWidth(0)[1];
		for (int i = 0; i < chain.length(); i++) {
			if ((this.row[i] - min_height) != (pr.row[i] - p_min_height)) {
				return false;
			}
			if ((this.col[i] - min_width) != (pr.col[i] - p_min_width)) {
				return false;
			}
		}
		return true;
	}

	/** Produces a random structure **/
	public void get_random_structure() {
		Protein p = (Protein) clone();
		// do some random number of global perturbations
		int num_global = ThreadLocalRandom.current().nextInt(1, this.chain.length() + 1);
		for (int i = 0; i < num_global; i++) {
			do {
				int rotations = ThreadLocalRandom.current().nextInt(1, 4);
				for (int j = 0; j < rotations; j++) {
					p.pivot_move();
				}
			} while (p.getScore() == -1);
			copy_arrays(p);
		}
		// do some random number of corner perterbations
		p.establish_corners();
		int total_corners = p.corners.size();
		if (total_corners == 0) {
			return;
		}
		int num_corners = ThreadLocalRandom.current().nextInt(1, total_corners + 1);
		for (int i = 0; i < num_corners; i++) {
			do {
				p.perterb_corner();
			} while (p.getScore() == -1);
			copy_arrays(p);
		}
		p.establish_flats();
		p.make_crankshaft();
		copy_arrays(p);
		// do some random number of crankshaft perterbations
		p.establish_crankshafts();
		int total_crankshafts = p.crankshafts.size();
		if (total_crankshafts == 0) {
			return;
		}
		int num_cranks = ThreadLocalRandom.current().nextInt(1, total_crankshafts + 1);
		for (int i = 0; i < num_cranks; i++) {
			do {
				p.perterb_crank();
			} while (p.getScore() == -1);
			copy_arrays(p);
		}
	}

	/** Walk function **/
	public void walk(int temp, int max_temp) {
		int global = ThreadLocalRandom.current().nextInt(0, max_temp);
		int corner = ThreadLocalRandom.current().nextInt(0, max_temp);
		int crank = ThreadLocalRandom.current().nextInt(0, max_temp);
		Protein p = (Protein) clone();
		if (global < (max_temp - temp)) {
			do {
				p.pivot_move();
			} while (p.getScore() == -1);
		}
		if (corner < temp) {
			p.establish_corners();
			do {
				p.perterb_corner();
			} while (p.getScore() == -1);
		}
		if (crank < (max_temp/4)) {
			p.establish_flats();
			p.make_crankshaft();
		}
		if (crank < (max_temp/4)) {
			p.establish_crankshafts();
			do {
				p.perterb_crank();
			} while (p.getScore() == -1);
		}
		copy_arrays(p);
	}

	public void copy_arrays(Protein p) {
		System.arraycopy(p.row, 0, this.row, 0, p.row.length);
		System.arraycopy(p.col, 0, this.col, 0, p.col.length);
	}

	/** Relocates the residue at index to specified y- and x-coordinate. */
	public void setLocation(int index, int rowLoc, int colLoc) {
		row[index] = rowLoc;
		col[index] = colLoc;
	}

	/** Creates a duplicate of this protein chain. */
	public Object clone() {
		try {
			Protein ret = (Protein) super.clone();
			ret.row = new int[this.row.length];
			System.arraycopy(this.row, 0, ret.row, 0, this.row.length);
			ret.col = new int[this.col.length];
			System.arraycopy(this.col, 0, ret.col, 0, this.col.length);
			return ret;
		} catch (CloneNotSupportedException e) {
			return new RuntimeException(e);
		}
	}

	/** Returns the score for this protein chain, which is
	 * negative if the chain overlaps itself and otherwise is
	 * the number of adjacent pairs of hydrophobic residues. */
	public int getScore() {
		int ret = 0;
		int length = chain.length();
		for (int i = 1; i < length; i++) {
			int dx = Math.abs(col[i] - col[i - 1]);
			int dy = Math.abs(row[i] - row[i - 1]);
			if (dx + dy != 1) {
				// chain is illegal since a chained pair is not adjacent on grid
				return -1;
			}
		}
		for (int i = 0; i < length; i++) {
			int r0 = row[i], c0 = col[i];
			for (int j = i + 1; j < length; j++) {
				int r1 = row[j], c1 = col[j];
				if (r0 == r1 && c0 == c1) {
					// chain is illegal since it wraps back onto itself
					return -1;
				}
				if (isHydrophobic(i) && isHydrophobic(j)) {
					if (r0 == r1) {
						if (c0 == c1 + 1 || c1 == c0 + 1) ret++;
					} else if (c0 == c1) {
						if (r0 == r1 + 1 || r1 == r0 + 1) ret++;
					}
				}
			}
		}
		return ret;
	}

	/** Displays the amino acid to the specified output stream. */
	public void print(PrintStream out) {
		// Compute the range of rows and columns represented.
		int minRow = Integer.MAX_VALUE;
		int maxRow = Integer.MIN_VALUE;
		int minCol = Integer.MAX_VALUE;
		int maxCol = Integer.MIN_VALUE;
		for (int i = 0; i < row.length; i++) {
			if (row[i] < minRow) minRow = row[i];
			if (row[i] > maxRow) maxRow = row[i];
			if (col[i] < minCol) minCol = col[i];
			if (col[i] > maxCol) maxCol = col[i];
		}

		// Create a grid representing what to display at each location.
		// (Each entry has 1 bit set if horizontal line connects to right
		// and 2 bit set if vertical line connects to below. If a residue
		// is at that location, its index is added once, multiplied by 4,
		// and placed into grid.)
		int[][] grid = new int[maxRow - minRow + 1][maxCol - minCol + 1];
		int r0 = row[0] - minRow;
		int c0 = col[0] - minCol;
		grid[r0][c0] = 4; // mark that amino 0 is found at (r0,c0)
		for (int i = 1; i < row.length; i++) {
			int r1 = row[i] - minRow;
			int c1 = col[i] - minCol;
			grid[r1][c1] += 4 * (i + 1); // mark than amino i is at (r1,c1)
			if (r0 == r1) {
				if (c1 == c0 + 1) grid[r0][c0] += 1;
				else if (c0 == c1 + 1) grid[r1][c1] += 1;
			} else {
				if (r1 == r0 + 1) grid[r0][c0] += 2;
				else if (r0 == r1 + 1) grid[r1][c1] += 2;
			}
			r0 = r1;
			c0 = c1;
		}

		// Display the grid as computed.
		for (int i = 0; i < grid.length; i++) {
			for (int j = 0; j < grid[0].length; j++) {
				int k = grid[i][j];
				int pos = k / 4 - 1;

				if (pos < 0) out.print(" ");
				else if (pos == 0) out.print(Character.toLowerCase(chain.charAt(pos)));
				else if (pos < chain.length()) out.print(chain.charAt(pos));
				else out.print("?");

				if (j < grid[0].length - 1) {
					if (k % 2 == 1) out.print("--");
					else out.print("  ");
				}
			}
			out.println();
			if (i < grid.length - 1) {
				for (int j = 0; j < grid[0].length; j++) {
					if (grid[i][j] % 4 >= 2) out.print("|");
					else out.print(" ");
					if (j < grid[0].length - 1) out.print("  ");
				}
				out.println();
			}
		}
	}

	public static void main(String[] args) {
		if (args.length == 0) {
			System.out.println("usage: java Protein CHAIN1 [CHAIN2 ...]");
			return;
		}
		long totalStart = System.nanoTime();
		int totalScore = 0;
		for (int arg = 0; arg < args.length; arg++) {
			best_score = 0;
			for (int i = 0; i < 10; i++) {
				totalScore += testChain(args[arg]);
			}
		}
		double totalElapse = (System.nanoTime() - totalStart) / 1e6;
		System.out.printf("%7.1f ms%3d TOTAL\n", totalElapse, totalScore);
	}

	private static int testChain(String chain) {
		if (chain.length() <= 1) {
			System.out.printf("%7.1f ms%3d %s\n", 0.0, 0, chain);
			return 0;
		} else {
			long start = System.nanoTime();
			Searcher searcher = new Searcher(chain);
			Protein best = searcher.search();
			double elapse = (System.nanoTime() - start) / 1e6;
			int score = best.getScore();
			best.print(System.out);
			if (best.getScore() > best_score) {
				searcher.graph_data();
				best_score = best.getScore();
			}
			System.out.println("Score:" + score);
			System.out.printf("%7.1f ms%3d %s\n", elapse, score, chain);
			return score;
		}
	}
}
