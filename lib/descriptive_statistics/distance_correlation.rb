module DescriptiveStatistics

  module DistanceCorrelation

    MISSING_POLICIES = [:drop, :raise].freeze

    module_function

    def clean_pairs(collection, missing)
      unless MISSING_POLICIES.include?(missing)
        raise ArgumentError, "unknown missing policy #{missing.inspect}; expected one of #{MISSING_POLICIES.inspect}"
      end

      xs = []
      ys = []
      incomplete_indexes = []
      collection.to_a.each_with_index do |row, index|
        pair = row.respond_to?(:to_a) ? row.to_a : nil
        unless pair && pair.size == 2
          raise ArgumentError, "row #{index} is not a pair of vectors: #{row.inspect}"
        end
        x = vector_or_nil(pair[0])
        y = vector_or_nil(pair[1])
        if x.nil? || y.nil?
          incomplete_indexes << index
        else
          xs << x
          ys << y
        end
      end

      if incomplete_indexes.any? && missing == :raise
        raise ArgumentError, "incomplete rows at indexes #{incomplete_indexes.inspect}"
      end

      check_dimensions(xs, ys)
      [xs, ys]
    end

    def vector_or_nil(value)
      return [value.to_f] if value.is_a?(Numeric)
      return nil unless value.respond_to?(:to_a)

      vector = value.to_a
      return nil if vector.empty?
      return nil if vector.any? { |element| missing_element?(element) }
      unless vector.all? { |element| element.is_a?(Numeric) }
        raise ArgumentError, "non-numeric vector element in #{vector.inspect}"
      end
      vector.map(&:to_f)
    end

    def missing_element?(element)
      element.nil? || (element.is_a?(Float) && element.nan?)
    end

    def check_dimensions(xs, ys)
      return if xs.empty?

      x_dimension = xs.first.size
      y_dimension = ys.first.size
      inconsistent = xs.any? { |x| x.size != x_dimension } ||
                     ys.any? { |y| y.size != y_dimension }
      if inconsistent
        raise ArgumentError, "inconsistent vector dimensions across rows"
      end
    end

    def require_size!(xs, unbiased)
      if unbiased
        if xs.size < 4
          raise ArgumentError, "unbiased distance covariance requires at least 4 complete pairs, got #{xs.size}"
        end
      elsif xs.size < 2
        raise ArgumentError, "distance covariance requires at least 2 complete pairs, got #{xs.size}"
      end
    end

    def distance_matrix(vectors)
      n = vectors.size
      matrix = Array.new(n) { Array.new(n, 0.0) }
      (0...n).each do |i|
        ((i + 1)...n).each do |j|
          distance = euclidean(vectors[i], vectors[j])
          matrix[i][j] = distance
          matrix[j][i] = distance
        end
      end
      matrix
    end

    def euclidean(u, v)
      sum = 0.0
      u.each_with_index do |component, k|
        difference = component - v[k]
        sum += difference * difference
      end
      Math.sqrt(sum)
    end

    def dcov_squared(xs, ys, unbiased)
      dcov_squared_from(distance_matrix(xs), distance_matrix(ys), unbiased)
    end

    def dcov_squared_from(a, b, unbiased)
      unbiased ? dcov_squared_u(a, b) : dcov_squared_v(a, b)
    end

    def dcov_squared_v(a, b)
      n = a.size
      a_row_means = a.map { |row| row.reduce(:+) / n }
      b_row_means = b.map { |row| row.reduce(:+) / n }
      a_grand_mean = a_row_means.reduce(:+) / n
      b_grand_mean = b_row_means.reduce(:+) / n
      total = 0.0
      n.times do |i|
        n.times do |j|
          total += (a[i][j] - a_row_means[i] - a_row_means[j] + a_grand_mean) *
                   (b[i][j] - b_row_means[i] - b_row_means[j] + b_grand_mean)
        end
      end
      total / (n * n)
    end

    def dcov_squared_u(a, b)
      n = a.size
      centered_a = u_center(a)
      centered_b = u_center(b)
      total = 0.0
      n.times do |i|
        n.times do |j|
          total += centered_a[i][j] * centered_b[i][j]
        end
      end
      total / (n * (n - 3))
    end

    def u_center(matrix)
      n = matrix.size
      row_sums = matrix.map { |row| row.reduce(:+) }
      grand_sum = row_sums.reduce(:+)
      centered = Array.new(n) { Array.new(n, 0.0) }
      n.times do |i|
        n.times do |j|
          next if i == j
          centered[i][j] = matrix[i][j] -
                           row_sums[i] / (n - 2.0) -
                           row_sums[j] / (n - 2.0) +
                           grand_sum / ((n - 1.0) * (n - 2.0))
        end
      end
      centered
    end

    def signed_square_root(value)
      value < 0 ? -Math.sqrt(-value) : Math.sqrt(value)
    end

  end

  def distance_covariance(collection = self, unbiased: false, missing: :drop)
    DistanceCorrelation.signed_square_root(
      distance_covariance_squared(collection, unbiased: unbiased, missing: missing))
  end

  def distance_covariance_squared(collection = self, unbiased: false, missing: :drop)
    xs, ys = DistanceCorrelation.clean_pairs(collection, missing)
    DistanceCorrelation.require_size!(xs, unbiased)
    DistanceCorrelation.dcov_squared(xs, ys, unbiased)
  end

  def distance_correlation(collection = self, unbiased: false, missing: :drop)
    DistanceCorrelation.signed_square_root(
      distance_correlation_squared(collection, unbiased: unbiased, missing: missing))
  end

  def distance_correlation_squared(collection = self, unbiased: false, missing: :drop)
    xs, ys = DistanceCorrelation.clean_pairs(collection, missing)
    DistanceCorrelation.require_size!(xs, unbiased)
    covariance = DistanceCorrelation.dcov_squared(xs, ys, unbiased)
    variance_x = DistanceCorrelation.dcov_squared(xs, xs, unbiased)
    variance_y = DistanceCorrelation.dcov_squared(ys, ys, unbiased)
    return 0.0 unless variance_x > 0.0 && variance_y > 0.0
    covariance / Math.sqrt(variance_x * variance_y)
  end

  def distance_covariance_permutation_test(collection = self, permutations: 999, seed: nil, unbiased: false, missing: :drop)
    unless permutations.is_a?(Integer) && permutations > 0
      raise ArgumentError, "permutations must be a positive Integer"
    end
    xs, ys = DistanceCorrelation.clean_pairs(collection, missing)
    DistanceCorrelation.require_size!(xs, unbiased)

    a = DistanceCorrelation.distance_matrix(xs)
    b = DistanceCorrelation.distance_matrix(ys)
    observed = DistanceCorrelation.dcov_squared_from(a, b, unbiased)
    rng = seed.nil? ? Random.new : Random.new(seed)
    n = xs.size
    exceedances = 0
    permutations.times do
      order = (0...n).to_a.shuffle(random: rng)
      permuted_b = order.map { |i| order.map { |j| b[i][j] } }
      statistic = DistanceCorrelation.dcov_squared_from(a, permuted_b, unbiased)
      exceedances += 1 if statistic >= observed
    end

    {
      statistic: observed,
      p_value: (1 + exceedances).to_f / (permutations + 1),
      permutations: permutations,
      estimator: unbiased ? :unbiased : :biased
    }
  end

end
