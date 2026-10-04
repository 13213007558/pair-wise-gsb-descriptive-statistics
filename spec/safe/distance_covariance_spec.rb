require 'safe/spec_helper'

describe "DescriptiveStatistics distance statistics" do
  require 'descriptive_statistics/safe'

  let(:pairs) do
    [
      [[0, 0], 0],
      [[1, 0], 1],
      [[2, 0], 2],
      [[3, 0], 3]
    ]
  end

  it "calculates the ordinary sample estimate from explicit 4x4 distance matrices" do
    x_distances = [
      [0.0, 1.0, 2.0, 3.0],
      [1.0, 0.0, 1.0, 2.0],
      [2.0, 1.0, 0.0, 1.0],
      [3.0, 2.0, 1.0, 0.0]
    ]
    y_distances = x_distances
    x_row_means = [1.5, 1.0, 1.0, 1.5]
    y_row_means = x_row_means
    grand_mean = 1.25

    x_centered = x_distances.each_with_index.map do |row, i|
      row.each_with_index.map do |distance, j|
        distance - x_row_means[i] - y_row_means[j] + grand_mean
      end
    end
    y_centered = x_centered
    squared_covariance = x_centered.flatten.zip(y_centered.flatten).map { |x, y| x * y }.reduce(:+) / 16.0

    statistics = DescriptiveStatistics.distance_dependence_statistics(pairs)

    expect(statistics[:squared_distance_covariance]).to eql(squared_covariance)
    expect(statistics[:distance_covariance]).to eql(Math.sqrt(squared_covariance))
    expect(statistics[:squared_distance_correlation]).to eql(1.0)
    expect(statistics[:distance_correlation]).to eql(1.0)
    expect(statistics[:squared_distance_variance_x]).to eql(squared_covariance)
    expect(statistics[:squared_distance_variance_y]).to eql(squared_covariance)
  end

  it "calculates the U-centered estimate from explicit 4x4 distance matrices" do
    distances = [
      [0.0, 1.0, 2.0, 3.0],
      [1.0, 0.0, 1.0, 2.0],
      [2.0, 1.0, 0.0, 1.0],
      [3.0, 2.0, 1.0, 0.0]
    ]
    row_sums = [6.0, 4.0, 4.0, 6.0]
    total = 20.0
    u_centered = Array.new(4) do |i|
      Array.new(4) do |j|
        if i == j
          0.0
        else
          distances[i][j] - (row_sums[i] - distances[i][j]) / 2.0 -
            (row_sums[j] - distances[i][j]) / 2.0 + total / 6.0
        end
      end
    end
    expected = u_centered.flatten.map { |value| value * value }.reduce(:+) / 4.0

    statistics = DescriptiveStatistics.distance_dependence_statistics(pairs)

    expect(statistics[:u_centered_squared_distance_covariance]).to eql(expected)
    expect(statistics[:u_centered_squared_distance_variance_x]).to eql(expected)
    expect(statistics[:u_centered_squared_distance_variance_y]).to eql(expected)
    expect(statistics[:u_centered_squared_distance_correlation]).to eql(1.0)
  end

  it "returns a negative U-centered squared estimate without clipping it to zero" do
    negative_pairs = [
      [[2, 0], 2],
      [[0, 0], -1],
      [[2, 0], -1],
      [[0, 0], 2],
      [[0, 0], 1]
    ]

    statistics = DescriptiveStatistics.distance_dependence_statistics(negative_pairs)

    expect(statistics[:u_centered_squared_distance_covariance]).to eql(-0.26666666666666666)
    expect(statistics[:u_centered_squared_distance_correlation]).to be_within(1e-15).of(-0.03481553119113956)
    expect(DescriptiveStatistics.distance_covariance(negative_pairs, :estimator => :u_centered_squared)).to eql(-0.26666666666666666)
    expect(DescriptiveStatistics.distance_correlation(negative_pairs, :estimator => :u_centered_squared)).to be_within(1e-15).of(-0.03481553119113956)
  end

  it "uses a fixed random source for paired permutation p values" do
    fixed_random = Class.new do
      attr_reader :calls

      def initialize
        @calls = 0
      end

      def rand(maximum)
        @calls += 1
        (@calls - 1) % maximum
      end
    end.new

    result = DescriptiveStatistics.distance_covariance_test(
      pairs,
      :repetitions => 9,
      :random => fixed_random
    )

    expect(result[:estimator]).to eql(:u_centered_squared)
    expect(result[:statistic]).to eql(12.000000000000004)
    expect(result[:p_value]).to eql(0.1)
    expect(result[:repetitions]).to eql(9)
    expect(result[:number]).to eql(4.0)
    expect(fixed_random.calls).to eql(27)
  end

  it "uses Random.new(seed) for reproducible permutation p values" do
    result = DescriptiveStatistics.distance_covariance_test(
      pairs,
      :repetitions => 23,
      :seed => 1234
    )

    expect(result[:p_value]).to eql(0.125)
  end

  it "uses the ordinary statistic for a sample estimator permutation test" do
    result = DescriptiveStatistics.distance_covariance_test(
      pairs,
      :estimator => :sample,
      :repetitions => 9,
      :seed => 1234
    )

    expect(result[:estimator]).to eql(:sample)
    expect(result[:statistic]).to eql(0.9013878188659973)
    expect(result[:p_value]).to eql(0.2)
  end

  it "supports separate paired collections and a block" do
    rows = [
      { :x => [0, 0], :y => 0 },
      { :x => [1, 0], :y => 1 }
    ]

    expect(DescriptiveStatistics.distance_covariance(rows) { |row| [row[:x], row[:y]] }).to eq(0.5)
    expect(DescriptiveStatistics.distance_covariance([0, 1, 2, 3], [0, 1, 2, 3])).to eq(0.9013878188659973)
  end

  it "uses values from a hash of paired rows" do
    rows = {
      :a => [[0, 0], 0],
      :b => [[1, 0], 1],
      :c => [[2, 0], 2],
      :d => [[3, 0], 3]
    }

    expect(DescriptiveStatistics.distance_correlation(rows)).to eql(1.0)
  end

  it "rejects separate collections with different numbers of rows" do
    expect do
      DescriptiveStatistics.distance_covariance([0, 1, 2, 3], [0, 1, 2])
    end.to raise_error(ArgumentError, "X and Y must contain the same number of paired rows")
  end

  it "supports options when the method is called on a collection" do
    expect(pairs.extend(DescriptiveStatistics).distance_correlation(:estimator => :u_centered_squared)).to eql(1.0)
  end

  it "defaults to dropping rows containing missing values and counts the dropped rows" do
    rows = [
      [[0, 0], 0],
      [[1, 0], nil],
      [[Float::NAN, 0], 2],
      [[3, 0], 3]
    ]

    statistics = DescriptiveStatistics.distance_dependence_statistics(rows)

    expect(statistics[:number]).to eql(2.0)
    expect(statistics[:dropped_rows]).to eql(2.0)
    expect(statistics[:distance_correlation]).to eql(1.0)
  end

  it "can raise when a paired row contains missing values" do
    rows = [[[0, 0], 0], [[1, 0], nil]]

    expect do
      DescriptiveStatistics.distance_dependence_statistics(rows, :missing_rows => :raise)
    end.to raise_error(ArgumentError, "missing value in pair at row 1")
  end

  it "accepts complete_cases as an alias for complete rows" do
    rows = [[[0, 0], 0], [[1, 0], nil], [[2, 0], 2], [[3, 0], 3]]

    statistics = DescriptiveStatistics.distance_dependence_statistics(rows, :missing_rows => :complete_cases)

    expect(statistics[:number]).to eql(3.0)
    expect(statistics[:dropped_rows]).to eql(1.0)
  end

  it "supports vectors with different dimensions for X and Y" do
    rows = [
      [[0, 0], [0, 0, 0]],
      [[1, 0], [1, 0, 0]],
      [[2, 0], [2, 0, 0]],
      [[3, 0], [3, 0, 0]]
    ]

    statistics = DescriptiveStatistics.distance_dependence_statistics(rows)

    expect(statistics[:distance_correlation]).to eql(1.0)
    expect(statistics[:u_centered_squared_distance_correlation]).to eql(1.0)
  end

  it "rejects rows whose vectors have inconsistent dimensions" do
    rows = [[[0, 0], 0], [[1, 0, 1], 1], [[2, 0], 2], [[3, 0], 3]]

    expect do
      DescriptiveStatistics.distance_dependence_statistics(rows)
    end.to raise_error(ArgumentError, /X vectors must all have the same dimension/)
  end

  it "returns zero statistics and a permutation p value of one for constant X" do
    constants = [[[1, 1], 8], [[1, 1], 9], [[1, 1], 10], [[1, 1], 11]]

    statistics = DescriptiveStatistics.distance_dependence_statistics(constants)
    test = DescriptiveStatistics.distance_covariance_test(constants, :repetitions => 9, :seed => 3)

    expect(statistics[:distance_covariance]).to eql(0.0)
    expect(statistics[:distance_correlation]).to eql(0.0)
    expect(statistics[:u_centered_squared_distance_covariance]).to eql(0.0)
    expect(statistics[:u_centered_squared_distance_correlation]).to eql(0.0)
    expect(test[:statistic]).to eql(0.0)
    expect(test[:p_value]).to eql(1.0)
  end

  it "requires four complete rows for the U-centered estimate" do
    rows = [[0, 0], [1, 1], [2, 2]]

    expect do
      DescriptiveStatistics.distance_covariance(rows, :estimator => :u_centered_squared)
    end.to raise_error(ArgumentError, "U-centered squared distance statistics require at least 4 complete rows")
  end
end
