// From http://www.cburch.com/ccsc-nifty/12/bur/handout.html
import java.util.ArrayList;
import java.util.concurrent.ThreadLocalRandom;
import java.lang.Math;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.data.category.CategoryDataset;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.awt.*;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import org.jfree.chart.ChartUtils;
import java.io.IOException;


public class Searcher extends JFrame {
	private Protein cur;
	private Protein best;
	private int bestScore;
	public int valids;
	public int invalids;
	public ArrayList<Protein> valid_proteins;
	int dups;
	int corners_turned_invalid;
	int crankshafts_turned_invalid;
	int pivot_turned_invalid;
	int[] scores;
	int max_temp;

	public Searcher(String chain) {
		cur = new Protein(chain);
		best = (Protein) cur.clone();
		bestScore = Integer.MIN_VALUE;
		valids = 0;
		invalids = 0;
		dups = 0;
		corners_turned_invalid = 0;
		crankshafts_turned_invalid = 0;
		pivot_turned_invalid = 0;
		valid_proteins = new ArrayList<Protein>();
		max_temp = 50000;
		scores = new int[max_temp + 1];
	}

	public Protein search() {
		int temp = 0;
		double e = Math.exp(1);
		cur.get_random_structure();
		int temp_of_first_best = 0;
		scores[temp] = cur.getScore();
		best = (Protein) cur.clone();
		Protein prev;
		long startTime = System.currentTimeMillis();
		while (true) {
			prev = (Protein) cur.clone();
			cur.walk(temp, max_temp);
			int score = cur.getScore();
			int prev_score = scores[temp];
			if (score > bestScore) {
				bestScore = score;
				best = (Protein) cur.clone();
				temp++;
				temp_of_first_best = temp;
				scores[temp] = score;
			} else if (score > prev_score) {
				temp++;
				scores[temp] = score;
			} else if (score < prev_score) {
				double range = (e/(bestScore - score))/(temp + 3);
				int size = ((int) ((temp + 3)/range));
				int val = ThreadLocalRandom.current().nextInt(0, size);
				double cap = (((double) range) * size);
				temp++;
				if (val > (cap - 1)) {
					cur = (Protein) prev.clone();
				}
				scores[temp] = cur.getScore();
			}
			if (temp == max_temp || ((System.currentTimeMillis() - startTime)/30000) > 5) {
				break;
			}
		}
		System.out.println(temp_of_first_best);
		return best;
	}

	public void graph_data() {
		XYSeriesCollection dataset = new XYSeriesCollection();
		XYSeries series = new XYSeries(cur.chain);
		for (int i = 0; i < max_temp; i++) {
			series.add(i, scores[i]);
		}
		dataset.addSeries(series);

		JFreeChart chart = ChartFactory.createXYLineChart(
			"Scores per Generation",	// Title
			"Iterations",			// X-axis label
			"Score",			// Y-axis label
			dataset
		);

		ChartPanel chartPanel = new ChartPanel(chart);
		chartPanel.setPreferredSize(new Dimension(560, 370));
		setContentPane(chartPanel);

		try {
			ChartUtils.saveChartAsPNG(new File(cur.chain), chart, 800, 600);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void printArray(int[] array) {
		for (int i = 0; i < 10; i++) {
			System.out.print(" " + array[i]);
		}
		System.out.println();
	}

	private void printMatrix(ArrayList<ArrayList<Integer>> matrix) {
		for (int i = 0; i < matrix.size(); i++) {
			for (int j = 0; j < matrix.get(i).size(); j++) {
				System.out.print(" " + matrix.get(i).get(j));
			}
			System.out.println();
		}
	}


}
