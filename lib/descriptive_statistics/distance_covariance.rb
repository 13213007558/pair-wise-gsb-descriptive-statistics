module DescriptiveStatistics

  DISTANCE_ESTIMATORS = [:sample, :u_centered_squared].freeze
  DISTANCE_MISSING_ROW_STRATEGIES = [:complete_rows, :complete_cases, :raise].freeze

  def distance_covariance(collection = self, y = nil, options = nil, &block)
    if collection.is_a?(Hash) && y.nil? && self.is_a?(Enumerable)
      options = collection
      collection = self
    end
    if y.is_a?(Hash) && options.nil?
      options = y
      y = nil
    end
    options ||= {}

    pairs = DescriptiveStatistics.send(:distance_pairs, collection, y, &block)
    analysis = DistanceAnalysis.new(pairs, options)
    analysis.distance_covariance
  end

  def distance_correlation(collection = self, y = nil, options = nil, &block)
    if collection.is_a?(Hash) && y.nil? && self.is_a?(Enumerable)
      options = collection
      collection = self
    end
    if y.is_a?(Hash) && options.nil?
      options = y
      y = nil
    end
    options ||= {}

    pairs = DescriptiveStatistics.send(:distance_pairs, collection, y, &block)
    analysis = DistanceAnalysis.new(pairs, options)
    analysis.distance_correlation
  end

  def distance_dependence_statistics(collection = self, y = nil, options = nil, &block)
    if collection.is_a?(Hash) && y.nil? && self.is_a?(Enumerable)
      options = collection
      collection = self
    end
    if y.is_a?(Hash) && options.nil?
      options = y
      y = nil
    end
    options ||= {}

    pairs = DescriptiveStatistics.send(:distance_pairs, collection, y, &block)
    DistanceAnalysis.new(pairs, options).statistics
  end

  def distance_covariance_test(collection = self, y = nil, options = nil, &block)
    if collection.is_a?(Hash) && y.nil? && self.is_a?(Enumerable)
      options = collection
      collection = self
    end
    if y.is_a?(Hash) && options.nil?
      options = y
      y = nil
    end
    options = options ? options.dup : {}
    options[:estimator] ||= :u_centered_squared

    pairs = DescriptiveStatistics.send(:distance_pairs, collection, y, &block)
    analysis = DistanceAnalysis.new(pairs, options)
    analysis.permutation_test
  end
  alias_method :distance_covariance_permutation_test, :distance_covariance_test

  private

  def distance_pairs(collection, y, &block)
    if block_given?
      rows = collection.map { |row| yield row }
    elsif collection.is_a?(Hash)
      rows = collection.values
    else
      rows = collection.to_a
    end
    return rows if y.nil?
    y_values = y.to_a
    if rows.size != y_values.size
      raise ArgumentError, "X and Y must contain the same number of paired rows"
    end
    rows.zip(y_values)
  end
  module_function :distance_pairs
  private_class_method :distance_pairs

  class DistanceAnalysis

    def initialize(pairs, options = {})
      @estimator = options.fetch(:estimator, :sample)
      @missing_rows = options.fetch(:missing_rows, :complete_rows)
      @missing_rows = :complete_rows if @missing_rows == :complete_cases
      @repetitions = options.fetch(:repetitions, 999)
      @seed = options[:seed]
      @random = options[:random]

      unless DISTANCE_ESTIMATORS.include?(@estimator)
        raise ArgumentError, "estimator must be one of #{DISTANCE_ESTIMATORS.map(&:inspect).join(', ')}"
      end
      unless DISTANCE_MISSING_ROW_STRATEGIES.include?(@missing_rows)
        raise ArgumentError, "missing_rows must be one of #{DISTANCE_MISSING_ROW_STRATEGIES.map(&:inspect).join(', ')}"
      end

      @x_values, @y_values, @dropped_rows = prepare_pairs(pairs)
      @sample_size = @x_values.size
      validate_sample_size
    end

    def distance_covariance
      if @estimator == :u_centered_squared
        u_statistics[:covariance]
      else
        Math.sqrt(sample_statistics[:squared_covariance])
      end
    end

    def distance_correlation
      if @estimator == :u_centered_squared
        u_statistics[:squared_correlation]
      else
        Math.sqrt(sample_statistics[:squared_correlation])
      end
    end

    def statistics
      sample = sample_statistics
      unbiased = u_statistics
      {
        :number => @sample_size.to_f,
        :dropped_rows => @dropped_rows.to_f,
        :distance_covariance => Math.sqrt(sample[:squared_covariance]),
        :distance_correlation => Math.sqrt(sample[:squared_correlation]),
        :squared_distance_covariance => sample[:squared_covariance],
        :squared_distance_correlation => sample[:squared_correlation],
        :squared_distance_variance_x => sample[:variance_x],
        :squared_distance_variance_y => sample[:variance_y],
        :u_centered_squared_distance_covariance => unbiased[:covariance],
        :u_centered_squared_distance_correlation => unbiased[:squared_correlation],
        :u_centered_squared_distance_variance_x => unbiased[:variance_x],
        :u_centered_squared_distance_variance_y => unbiased[:variance_y]
      }
    end

    def permutation_test
      repetitions = Integer(@repetitions)
      raise ArgumentError, "repetitions must be non-negative" if repetitions < 0

      observed = permutation_statistics[:covariance]
      random = @random || Random.new(@seed)
      exceedances = 0

      repetitions.times do
        permuted_y = permute(@y_values, random)
        permuted = DistanceAnalysis.new(@x_values.zip(permuted_y), permutation_analysis_options)
        exceedances += 1 if permuted.permutation_statistics[:covariance] >= observed
      end

      {
        :estimator => @estimator,
        :statistic => observed,
        :p_value => (exceedances + 1).to_f / (repetitions + 1),
        :repetitions => repetitions,
        :number => @sample_size.to_f,
        :dropped_rows => @dropped_rows.to_f
      }
    end

    protected

    def permutation_statistics
      if @estimator == :u_centered_squared
        u_statistics
      else
        sample_statistics.merge(:covariance => Math.sqrt(sample_statistics[:squared_covariance]))
      end
    end

    private

    def prepare_pairs(pairs)
      rows = pairs.to_a
      complete_rows = []
      dropped_rows = 0

      rows.each_with_index do |row, index|
        validate_pair_shape(row, index)
        x_vector = prepare_vector(row[0], index, 0)
        y_vector = prepare_vector(row[1], index, 1)
        missing = vector_missing?(x_vector) || vector_missing?(y_vector)

        if missing
          dropped_rows += 1
          if @missing_rows == :raise
            raise ArgumentError, "missing value in pair at row #{index}"
          end
        else
          validate_vector_values(x_vector, index, 0)
          validate_vector_values(y_vector, index, 1)
          complete_rows << [x_vector, y_vector]
        end
      end

      x_values = complete_rows.map { |row| row[0] }
      y_values = complete_rows.map { |row| row[1] }
      validate_dimensions(x_values, 'X')
      validate_dimensions(y_values, 'Y')
      [x_values, y_values, dropped_rows]
    end

    def validate_pair_shape(row, index)
      unless row.is_a?(Array) && row.size == 2
        raise ArgumentError, "pair at row #{index} must be an Array containing two vectors"
      end
    end

    def prepare_vector(value, index, side)
      return [nil] if value.nil?
      return [value] if value.is_a?(Numeric)
      unless value.is_a?(Array)
        raise ArgumentError, "#{side == 0 ? 'X' : 'Y'} vector at row #{index} must be a Numeric or an Array of Numerics"
      end
      value
    end

    def vector_missing?(vector)
      vector.any? { |value| nil == value || (value.is_a?(Float) && value.nan?) }
    end

    def validate_vector_values(vector, index, side)
      vector.each do |value|
        unless value.is_a?(Numeric) && (!value.respond_to?(:infinite?) || value.infinite?.nil?)
          raise ArgumentError, "#{side == 0 ? 'X' : 'Y'} vector at row #{index} contains a non-finite or non-Numeric value"
        end
      end
      if vector.empty?
        raise ArgumentError, "#{side == 0 ? 'X' : 'Y'} vector at row #{index} is empty"
      end
    end

    def validate_dimensions(values, name)
      return if values.empty?
      dimension = values.first.size
      values.each_with_index do |vector, index|
        unless vector.size == dimension
          raise ArgumentError, "#{name} vectors must all have the same dimension; row #{index} has #{vector.size} where #{dimension} was expected"
        end
      end
    end

    def validate_sample_size
      minimum_size = @estimator == :u_centered_squared ? 4 : 2
      if @x_values.size < minimum_size
        raise ArgumentError, "#{@estimator == :u_centered_squared ? 'U-centered squared' : 'sample'} distance statistics require at least #{minimum_size} complete rows"
      end
    end

    def sample_statistics
      @sample_statistics ||= centered_statistics(sample_centered(distance_matrix(@x_values)), sample_centered(distance_matrix(@y_values)), @sample_size * @sample_size.to_f)
    end

    def u_statistics
      @u_statistics ||= begin
        x_u = u_centered(distance_matrix(@x_values))
        y_u = u_centered(distance_matrix(@y_values))
        statistics = centered_statistics(x_u, y_u, @sample_size * (@sample_size - 3).to_f)
        statistics[:squared_correlation] = u_squared_correlation(statistics)
        statistics
      end
    end

    def centered_statistics(x_centered, y_centered, divisor)
      covariance = inner_product(x_centered, y_centered) / divisor
      variance_x = inner_product(x_centered, x_centered) / divisor
      variance_y = inner_product(y_centered, y_centered) / divisor
      denominator = Math.sqrt(variance_x * variance_y)
      squared_correlation = denominator.zero? ? 0.0 : covariance / denominator
      {
        :covariance => covariance,
        :squared_covariance => covariance,
        :variance_x => variance_x,
        :variance_y => variance_y,
        :squared_correlation => squared_correlation
      }
    end

    def u_squared_correlation(statistics)
      denominator = Math.sqrt(statistics[:variance_x] * statistics[:variance_y])
      return 0.0 if denominator.zero?
      statistics[:covariance] / denominator
    end

    def distance_matrix(values)
      values.map do |left|
        values.map do |right|
          Math.sqrt(left.zip(right).map { |a, b| (a - b) ** 2 }.reduce(0.0, :+))
        end
      end
    end

    def sample_centered(matrix)
      row_means = matrix.map { |row| row.reduce(0.0, :+) / @sample_size }
      grand_mean = row_means.reduce(0.0, :+) / @sample_size

      matrix.each_with_index.map do |row, i|
        row.each_with_index.map do |distance, j|
          distance - row_means[i] - row_means[j] + grand_mean
        end
      end
    end

    def u_centered(matrix)
      row_sums = matrix.map { |row| row.reduce(0.0, :+) }
      total = row_sums.reduce(0.0, :+)

      Array.new(@sample_size) do |i|
        Array.new(@sample_size) do |j|
          if i == j
            0.0
          else
            matrix[i][j] -
              (row_sums[i] - matrix[i][j]) / (@sample_size - 2.0) -
              (row_sums[j] - matrix[i][j]) / (@sample_size - 2.0) +
              total / ((@sample_size - 1.0) * (@sample_size - 2.0))
          end
        end
      end
    end

    def inner_product(left, right)
      left.each_with_index.reduce(0.0) do |total, (row, i)|
        total + row.each_with_index.reduce(0.0) do |row_total, (_, j)|
          row_total + left[i][j] * right[i][j]
        end
      end
    end

    def permute(values, random)
      permuted = values.dup
      (permuted.size - 1).downto(1) do |i|
        j = random.rand(i + 1)
        permuted[i], permuted[j] = permuted[j], permuted[i]
      end
      permuted
    end

    def permutation_analysis_options
      {
        :estimator => @estimator,
        :missing_rows => @missing_rows
      }
    end

  end

  private_constant :DistanceAnalysis

end
