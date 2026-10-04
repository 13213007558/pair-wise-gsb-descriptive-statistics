require 'monkeypatch/spec_helper'

describe "distance covariance and correlation" do
  require 'descriptive_statistics'

  let(:pairs3) { [[[0.0], [0.0]], [[1.0], [2.0]], [[3.0], [1.0]]] }
  let(:pairs4) { [[[0.0], [0.0]], [[1.0], [1.0]], [[2.0], [2.0]], [[4.0], [4.0]]] }

  def expanded_distance_matrix(vectors)
    vectors.map do |u|
      vectors.map do |v|
        Math.sqrt(u.zip(v).map { |a, b| (a - b)**2 }.reduce(0.0, :+))
      end
    end
  end

  def expanded_dcov_squared_v(xs, ys)
    a = expanded_distance_matrix(xs)
    b = expanded_distance_matrix(ys)
    n = xs.size
    a_row = a.map { |row| row.reduce(:+) / n }
    b_row = b.map { |row| row.reduce(:+) / n }
    a_grand = a_row.reduce(:+) / n
    b_grand = b_row.reduce(:+) / n
    total = 0.0
    n.times do |i|
      n.times do |j|
        total += (a[i][j] - a_row[i] - a_row[j] + a_grand) *
                 (b[i][j] - b_row[i] - b_row[j] + b_grand)
      end
    end
    total / (n * n)
  end

  def expanded_u_center(matrix)
    n = matrix.size
    row_sums = matrix.map { |row| row.reduce(:+) }
    grand = row_sums.reduce(:+)
    Array.new(n) do |i|
      Array.new(n) do |j|
        if i == j
          0.0
        else
          matrix[i][j] - row_sums[i] / (n - 2.0) - row_sums[j] / (n - 2.0) + grand / ((n - 1.0) * (n - 2.0))
        end
      end
    end
  end

  def expanded_dcov_squared_u(xs, ys)
    a = expanded_u_center(expanded_distance_matrix(xs))
    b = expanded_u_center(expanded_distance_matrix(ys))
    n = xs.size
    total = 0.0
    n.times { |i| n.times { |j| total += a[i][j] * b[i][j] } }
    total / (n * (n - 3))
  end

  context "biased (V-statistic) estimates checked against hand-expanded distance matrices" do
    it "matches the hand-computed squared distance covariance for three rows" do
      expect(pairs3.distance_covariance_squared).to be_within(1e-12).of(4.0 / 9)
    end

    it "matches the hand-computed distance covariance" do
      expect(pairs3.distance_covariance).to be_within(1e-12).of(2.0 / 3)
    end

    it "matches the hand-computed squared distance correlation" do
      expect(pairs3.distance_correlation_squared).to be_within(1e-12).of(9.0 / (4 * Math.sqrt(15)))
    end

    it "matches the hand-computed distance correlation" do
      expect(pairs3.distance_correlation).to be_within(1e-12).of(Math.sqrt(9.0 / (4 * Math.sqrt(15))))
    end

    it "agrees with an independent matrix expansion on multivariate rows of different dimensions" do
      pairs = [
        [[0.0, 1.0], [0.0]],
        [[1.0, 0.0], [2.0]],
        [[3.0, 3.0], [1.0]],
        [[2.0, 1.0], [5.0]],
        [[0.0, 0.0], [3.0]]
      ]
      xs = pairs.map(&:first)
      ys = pairs.map(&:last)
      expect(pairs.distance_covariance_squared).to be_within(1e-12).of(expanded_dcov_squared_v(xs, ys))
      expect(pairs.distance_covariance_squared(unbiased: true)).to be_within(1e-12).of(expanded_dcov_squared_u(xs, ys))
    end
  end

  context "unbiased (U-centered) squared estimates" do
    it "matches the hand-computed unbiased squared distance covariance for four rows" do
      expect(pairs4.distance_covariance_squared(unbiased: true)).to be_within(1e-12).of(2.0 / 3)
    end

    it "matches the hand-computed unbiased squared distance correlation for four rows" do
      expect(pairs4.distance_correlation_squared(unbiased: true)).to be_within(1e-12).of(1.0)
    end

    it "returns a negative unbiased squared estimate as-is instead of clipping it to zero" do
      alternating = [[[0.0], [0.0]], [[0.0], [1.0]], [[1.0], [0.0]], [[1.0], [1.0]]]
      expect(alternating.distance_covariance_squared(unbiased: true)).to be_within(1e-12).of(-1.0 / 3)
      expect(alternating.distance_correlation_squared(unbiased: true)).to be_within(1e-12).of(-0.5)
    end

    it "reports the sign of a negative squared estimate through a signed square root" do
      alternating = [[[0.0], [0.0]], [[0.0], [1.0]], [[1.0], [0.0]], [[1.0], [1.0]]]
      expect(alternating.distance_covariance(unbiased: true)).to be_within(1e-12).of(-Math.sqrt(1.0 / 3))
      expect(alternating.distance_correlation(unbiased: true)).to be_within(1e-12).of(-Math.sqrt(0.5))
    end

    it "differs from the biased estimate on the same data" do
      expect(pairs4.distance_covariance_squared(unbiased: true)).not_to be_within(1e-12).of(pairs4.distance_covariance_squared)
    end

    it "requires at least four complete pairs" do
      expect { pairs3.distance_covariance_squared(unbiased: true) }.to raise_error(ArgumentError)
    end
  end

  context "detecting nonlinear dependence where Pearson is zero" do
    let(:quadratic) { (-10..10).map { |i| x = i / 10.0; [[x], [x * x]] } }

    it "finds a clearly positive distance correlation" do
      expect(quadratic.distance_correlation).to be > 0.2
    end

    it "has near-zero Pearson correlation on the same rows" do
      xs = quadratic.map { |row| row[0][0] }
      ys = quadratic.map { |row| row[1][0] }
      mean_x = xs.reduce(:+) / xs.size
      mean_y = ys.reduce(:+) / ys.size
      covariance = xs.zip(ys).map { |x, y| (x - mean_x) * (y - mean_y) }.reduce(:+)
      pearson = covariance / Math.sqrt(xs.map { |x| (x - mean_x)**2 }.reduce(:+) * ys.map { |y| (y - mean_y)**2 }.reduce(:+))
      expect(pearson.abs).to be < 1e-12
    end
  end

  context "missing rows" do
    let(:with_missing) do
      pairs3 + [
        [nil, [5.0]],
        [[4.0], nil],
        [[Float::NAN, 1.0], [1.0]],
        [[], [1.0]],
        [[2.0], []]
      ]
    end

    it "drops incomplete rows by default" do
      expect(with_missing.distance_covariance).to eql(pairs3.distance_covariance)
      expect(with_missing.distance_correlation).to eql(pairs3.distance_correlation)
    end

    it "raises with missing: :raise" do
      expect { with_missing.distance_covariance(missing: :raise) }.to raise_error(ArgumentError, /incomplete rows/)
    end

    it "rejects an unknown missing policy" do
      expect { pairs3.distance_covariance(missing: :ignore) }.to raise_error(ArgumentError, /missing policy/)
    end
  end

  context "input validation" do
    it "accepts scalars as one-dimensional vectors" do
      expect([[0, 0], [1, 2], [3, 1]].distance_covariance_squared).to be_within(1e-12).of(4.0 / 9)
    end

    it "rejects inconsistent vector dimensions across rows" do
      expect { [[[0.0, 1.0], [0.0]], [[1.0], [1.0]]].distance_correlation }.to raise_error(ArgumentError, /dimensions/)
    end

    it "rejects non-numeric vector elements" do
      expect { [[["a"], [0.0]], [[1.0], [1.0]]].distance_correlation }.to raise_error(ArgumentError, /non-numeric/)
    end

    it "rejects rows that are not pairs" do
      expect { [[[0.0], [0.0], [1.0]], [[1.0], [1.0]]].distance_correlation }.to raise_error(ArgumentError, /pair/)
    end

    it "requires at least two complete pairs" do
      expect { [[[0.0], [0.0]]].distance_covariance }.to raise_error(ArgumentError, /at least 2/)
      expect { [].distance_covariance }.to raise_error(ArgumentError, /at least 2/)
    end
  end

  context "degenerate constant input" do
    let(:constant_x) { [[[1.0], [0.0]], [[1.0], [1.0]], [[1.0], [2.0]], [[1.0], [3.0]]] }
    let(:constant_both) { [[[1.0], [2.0]], [[1.0], [2.0]], [[1.0], [2.0]], [[1.0], [2.0]]] }

    it "has zero distance covariance and correlation" do
      expect(constant_x.distance_covariance).to eql(0.0)
      expect(constant_x.distance_correlation).to eql(0.0)
      expect(constant_x.distance_covariance(unbiased: true)).to eql(0.0)
      expect(constant_x.distance_correlation(unbiased: true)).to eql(0.0)
      expect(constant_both.distance_correlation).to eql(0.0)
    end

    it "has a permutation p-value of one" do
      result = constant_x.distance_covariance_permutation_test(permutations: 99, seed: 42)
      expect(result[:statistic]).to eql(0.0)
      expect(result[:p_value]).to eql(1.0)
    end
  end

  context "permutation test" do
    let(:dependent) { (1..30).map { |i| [[i.to_f], [2.0 * i]] } }
    let(:independentish) { (1..30).map { |i| [[i.to_f], [(i * 7 % 13).to_f]] } }

    it "is reproducible with a fixed seed" do
      first = dependent.distance_covariance_permutation_test(permutations: 199, seed: 1234)
      second = dependent.distance_covariance_permutation_test(permutations: 199, seed: 1234)
      expect(first).to eql(second)
    end

    it "reports the estimator and permutation count" do
      result = dependent.distance_covariance_permutation_test(permutations: 99, seed: 5)
      expect(result[:permutations]).to eql(99)
      expect(result[:estimator]).to eql(:biased)
      expect(dependent.distance_covariance_permutation_test(permutations: 99, seed: 5, unbiased: true)[:estimator]).to eql(:unbiased)
    end

    it "uses the squared distance covariance of the observed pairing as the statistic" do
      result = dependent.distance_covariance_permutation_test(permutations: 99, seed: 5)
      expect(result[:statistic]).to eql(dependent.distance_covariance_squared)
    end

    it "is small for strongly dependent data" do
      result = dependent.distance_covariance_permutation_test(permutations: 199, seed: 1234)
      expect(result[:p_value]).to be <= 0.05
    end

    it "is not significant for unrelated data" do
      result = independentish.distance_covariance_permutation_test(permutations: 199, seed: 99)
      expect(result[:p_value]).to be > 0.05
    end

    it "works with the unbiased estimator" do
      result = dependent.distance_covariance_permutation_test(permutations: 199, seed: 1234, unbiased: true)
      expect(result[:statistic]).to eql(dependent.distance_covariance_squared(unbiased: true))
      expect(result[:p_value]).to be <= 0.05
    end

    it "rejects a non-positive permutation count" do
      expect { dependent.distance_covariance_permutation_test(permutations: 0) }.to raise_error(ArgumentError, /permutations/)
    end
  end

  context "module function access" do
    it "is available on DescriptiveStatistics" do
      expect(DescriptiveStatistics.distance_correlation(pairs3)).to be_within(1e-12).of(pairs3.distance_correlation)
      expect(DescriptiveStatistics.distance_covariance_squared(pairs4, unbiased: true)).to be_within(1e-12).of(2.0 / 3)
    end
  end
end
