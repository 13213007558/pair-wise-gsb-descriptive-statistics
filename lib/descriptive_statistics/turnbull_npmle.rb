module DescriptiveStatistics
  class TurnbullNpmle
    TYPES = [:exact, :left, :right, :interval].freeze

    Result = Struct.new(:support, :probabilities, :iterations, :log_likelihood,
                        :converged, :tolerance, :support_size, :active_set_size,
                        :active_set_iterations, :kkt_residual) do
      def cdf_range(point)
        point = number(point)
        lower = 0.0
        upper = 0.0

        support.each_with_index do |region, index|
          probability = probabilities[index].to_f
          upper += probability if point_strictly_above_region?(region, point)
          lower += probability if region_weakly_below_point?(region, point)
        end

        lower = 1.0 if lower > 1.0
        upper = 1.0 if upper > 1.0
        {lower: lower, upper: upper}
      end

      def cdf(point)
        range = cdf_range(point)
        range[:lower] if (range[:lower] - range[:upper]).abs <= 1e-15
      end

      private

      def number(value)
        Float(value)
      rescue ArgumentError, TypeError
        raise ArgumentError, "CDF point must be numeric"
      end

      def point_strictly_above_region?(region, point)
        lower = region[:lower]
        return true if lower.infinite?
        region[:left_closed] ? point >= lower : point > lower
      end

      def region_weakly_below_point?(region, point)
        upper = region[:upper]
        return false if upper.infinite?
        point >= upper
      end
    end

    SupportRegion = Struct.new(:lower, :upper, :left_closed, :right_closed,
                               :identifiable) do
      def to_h
        {
          lower: lower,
          upper: upper,
          left_closed: left_closed,
          right_closed: right_closed,
          identifiable: identifiable
        }
      end
    end

    def initialize(tolerance: 1e-10, max_iterations: 10_000,
                   zero_tolerance: 1e-14, kkt_tolerance: 1e-8)
      @tolerance = Float(tolerance)
      @max_iterations = Integer(max_iterations)
      @zero_tolerance = Float(zero_tolerance)
      @kkt_tolerance = Float(kkt_tolerance)
    end

    def fit(records)
      observations = records.map { |record| normalize(record) }
      raise ArgumentError, "at least one observation is required" if observations.empty?

      groups = build_support_groups(observations)
      raise ArgumentError, "observations do not cover any common support region" if groups.empty?

      indicators = build_indicators(observations, groups)
      probabilities, iterations, log_likelihood, converged =
        expectation_maximization(indicators, groups.size, 0...groups.size)

      active = probabilities.each_index.to_a
      active_iterations = 0
      active_converged = false

      groups.size.times do
        break if active.empty?
        scores = column_scores(indicators, probabilities)
        positive = active.select { |index| probabilities[index] > 1e-8 }
        target = positive.empty? ? nil : scores.values_at(*positive).sum / positive.size.to_f
        removable = active.select do |index|
          probabilities[index] < 1e-8 &&
            (positive.empty? || scores[index] <= target + @kkt_tolerance)
        end
        if removable.empty?
          active_converged = converged
          break
        end
        active -= removable unless removable.empty?

        active_probabilities, run_iterations, log_likelihood, converged =
          constrained_newton(indicators, groups.size, active)
        active_iterations += run_iterations
        active_converged = converged
        probabilities = active_probabilities
      end

      normalize_probabilities!(probabilities)
      active_indices = probabilities.each_index.select do |index|
        probabilities[index] > @zero_tolerance
      end
      log_likelihood = log_likelihood(indicators, probabilities)
      kkt_residual = kkt_residual(indicators, probabilities, active_indices)

      positive_indices = active_indices
      support = positive_indices.map do |index|
        group = groups[index]
        first = group.first
        last = group.last
        SupportRegion.new(
          first[:lower],
          last[:upper],
          first[:left_closed],
          last[:right_closed],
          group.size == 1 && first[:type] == :singleton
        )
      end

      Result.new(
        support.map(&:to_h),
        positive_indices.map { |index| probabilities[index] },
        iterations,
        log_likelihood,
        active_converged && kkt_residual <= @kkt_tolerance,
        @tolerance,
        groups.size,
        active.size,
        active_iterations,
        kkt_residual
      )
    end

    private

    def normalize(record)
      raise ArgumentError, "observation must be a Hash" unless record.is_a?(Hash)

      type = record[:type]
      lower = nil
      upper = nil
      left_closed = nil
      right_closed = nil

      case type
      when :exact
        value = finite_number(record[:value])
        lower = upper = value
        left_closed = right_closed = true
      when :left
        value = finite_number(record[:value])
        lower = -Float::INFINITY
        upper = value
        left_closed = false
        right_closed = required_bool(record, :right_closed)
      when :right
        value = finite_number(record[:value])
        lower = value
        upper = Float::INFINITY
        left_closed = required_bool(record, :left_closed)
        right_closed = false
      when :interval
        lower = number(record[:lower])
        upper = number(record[:upper])
        raise ArgumentError, "interval lower endpoint must be below upper endpoint" if lower >= upper
        left_closed = required_bool(record, :left_closed)
        right_closed = required_bool(record, :right_closed)
        validate_infinite_endpoints(lower, upper, left_closed, right_closed)
      when nil
        if record.key?(:value)
          value = finite_number(record[:value])
          lower = upper = value
          left_closed = right_closed = true
          type = :exact
        else
          lower = number(record[:lower])
          upper = number(record[:upper])
          raise ArgumentError, "interval lower endpoint must be below upper endpoint" if lower >= upper
          left_closed = required_bool(record, :left_closed)
          right_closed = required_bool(record, :right_closed)
          validate_infinite_endpoints(lower, upper, left_closed, right_closed)
          type = :interval
        end
      else
        raise ArgumentError, "unknown observation type: #{type.inspect}"
      end

      {type: type, lower: lower, upper: upper,
       left_closed: left_closed, right_closed: right_closed}
    end

    def bool(value, name)
      case value
      when true then true
      when false then false
      else raise ArgumentError, "#{name} must be true or false"
      end
    end

    def required_bool(record, name)
      unless record.key?(name)
        raise ArgumentError, "#{name} must be declared explicitly"
      end
      bool(record[name], name)
    end

    def validate_infinite_endpoints(lower, upper, left_closed, right_closed)
      if lower.infinite? && left_closed
        raise ArgumentError, "left_closed cannot be true at an infinite endpoint"
      end
      if upper.infinite? && right_closed
        raise ArgumentError, "right_closed cannot be true at an infinite endpoint"
      end
    end

    def number(value)
      Float(value)
    rescue ArgumentError, TypeError, NoMethodError
      raise ArgumentError, "endpoint must be numeric"
    end

    def finite_number(value)
      result = number(value)
      raise ArgumentError, "exact and censoring values must be finite" unless result.finite?
      result
    end

    def build_support_groups(observations)
      points = observations.flat_map do |observation|
        [observation[:lower], observation[:upper]].select(&:finite?)
      end.uniq.sort

      cells = []
      add_cell(cells, point: nil, type: :ray_left,
               lower: -Float::INFINITY, upper: points.first,
               left_closed: false, right_closed: false) unless points.empty?
      add_cell(cells, point: nil, type: :ray_right,
               lower: points.last, upper: Float::INFINITY,
               left_closed: false, right_closed: false) unless points.empty?

      points.each do |point|
        add_cell(cells, point: point, type: :singleton, lower: point, upper: point,
                 left_closed: true, right_closed: true)
      end

      points.each_cons(2) do |left_point, right_point|
        add_cell(cells, point: nil, type: :interval, lower: left_point,
                 upper: right_point, left_closed: false, right_closed: false)
      end

      if points.empty?
        add_cell(cells, point: nil, type: :interval,
                 lower: -Float::INFINITY, upper: Float::INFINITY,
                 left_closed: false, right_closed: false)
      end

      feasible = cells.sort { |left, right| cell_order(left, right) }.map do |cell|
        mask = observations.each_with_index.inject(0) do |mask, (observation, index)|
          cell_in_observation?(cell, observation) ? mask | (1 << index) : mask
        end
        next nil if mask.zero?
        cell.merge(mask: mask)
      end

      groups = []
      current = []

      feasible.each do |cell|
        if cell.nil?
          groups << current unless current.empty?
          current = []
        elsif !current.empty? && cell[:mask] == current.last[:mask] && extend?(current.last, cell)
          current << cell
        else
          groups << current unless current.empty?
          current = [cell]
        end
      end
      groups << current unless current.empty?
      groups
    end

    def add_cell(cells, **cell)
      cells << cell
    end

    def cell_order(left, right)
      left_key = left[:lower].infinite? ? -Float::INFINITY : left[:lower]
      right_key = right[:lower].infinite? ? -Float::INFINITY : right[:lower]
      return left_key <=> right_key unless left_key == right_key
      type_rank(left[:type]) <=> type_rank(right[:type])
    end

    def type_rank(type)
      {singleton: 0, interval: 1, ray_left: 1, ray_right: 1}.fetch(type)
    end

    def extend?(previous, cell)
      return false unless previous
      return false unless previous[:upper] == cell[:lower]
      previous[:right_closed] || cell[:left_closed]
    end

    def cell_in_observation?(cell, observation)
      if cell[:type] == :singleton
        point = cell[:lower]
        return point_in_observation?(point, observation)
      end

      lower = cell[:lower]
      upper = cell[:upper]
      lower >= observation[:lower] && upper <= observation[:upper]
    end

    def point_in_observation?(point, observation)
      lower = observation[:lower]
      upper = observation[:upper]
      return false unless point > lower || point == lower && observation[:left_closed]
      return false unless point < upper || point == upper && observation[:right_closed]
      true
    end

    def build_indicators(observations, groups)
      observations.map do |observation|
        groups.map do |group|
          group.any? { |cell| cell_in_observation?(cell, observation) }
        end
      end
    end

    def expectation_maximization(indicators, size, active_indices)
      probabilities = Array.new(size, 0.0)
      active_indices.each { |index| probabilities[index] = 1.0 / active_indices.size }
      previous_log = log_likelihood(indicators, probabilities)
      converged = false
      iterations = 0

      @max_iterations.times do
        iterations += 1
        denominators = Array.new(size, 0.0)
        indicators.each do |row|
          row_probability = 0.0
          row.each_with_index do |included, index|
            row_probability += probabilities[index] if included
          end
          next if row_probability.zero?

          row.each_with_index do |included, index|
            denominators[index] += probabilities[index] / row_probability if included
          end
        end

        updated = Array.new(size, 0.0)
        active_indices.each do |index|
          updated[index] = denominators[index] / indicators.size.to_f
        end
        normalize_probabilities!(updated)
        current_log = log_likelihood(indicators, updated)
        max_change = updated.each_with_index.map do |probability, index|
          (probability - probabilities[index]).abs
        end.max
        probabilities = updated
        previous_log = current_log
        if max_change <= @tolerance
          converged = true
          break
        end
      end

      [probabilities, iterations, previous_log, converged]
    end

    def log_likelihood(indicators, probabilities)
      indicators.sum do |row|
        row_probability = 0.0
        row.each_with_index do |included, index|
          row_probability += probabilities[index] if included
        end
        raise ArgumentError, "no support region covers an observation" if row_probability.zero?
        Math.log(row_probability)
      end
    end

    def constrained_newton(indicators, size, active_indices)
      probabilities = Array.new(size, 0.0)
      active_indices.each { |index| probabilities[index] = 1.0 / active_indices.size }
      previous_log = log_likelihood(indicators, probabilities)
      converged = false
      iterations = 0

      @max_iterations.times do
        iterations += 1
        gradient = column_scores(indicators, probabilities)
        hessian = active_indices.each_with_index.map do |index, row_index|
          active_indices.each_with_index.map do |other, column_index|
            value = indicators.inject(0.0) do |total, row|
              if row[index] && row[other]
                total - 1.0 / (row_probability(row, probabilities)**2)
              else
                total
              end
            end
            value
          end
        end
        rhs = active_indices.map { |index| -gradient[index] } + [0.0]
        kkt = hessian.map { |row| row + [1.0] }
        kkt << active_indices.map { 1.0 } + [0.0]
        solution = begin
          solve_symmetric(kkt, rhs)
        rescue ArgumentError
          regularized_hessian = hessian.each_with_index.map do |row, row_index|
            row.each_with_index.map do |value, column_index|
              value + (row_index == column_index ? 1e-12 : 0.0)
            end
          end
          regularized_kkt = regularized_hessian.map { |row| row + [1.0] }
          regularized_kkt << active_indices.map { 1.0 } + [0.0]
          solve_symmetric(regularized_kkt, rhs)
        end
        step = solution[0...active_indices.size]
        max_negative = step.each_with_index.map do |value, offset|
          value < 0.0 ? -value / probabilities[active_indices[offset]] : 0.0
        end.max
        step_size = max_negative.nil? || max_negative.zero? ? 1.0 : 0.99 / max_negative
        candidate = nil
        80.times do
          candidate = Array.new(size, 0.0)
          step.each_with_index do |value, offset|
            candidate[active_indices[offset]] =
              probabilities[active_indices[offset]] + step_size * value
          end
          feasible = candidate.values_at(*active_indices).all? { |value| value >= 0.0 }
          directional_derivative = step.each_with_index.inject(0.0) do |total, (value, offset)|
            total + value * gradient[active_indices[offset]]
          end
          expected_increase = 1e-4 * step_size * directional_derivative
          if feasible &&
             log_likelihood(indicators, candidate) >= previous_log + expected_increase
            probabilities = candidate
            break
          end
          step_size *= 0.5
        end
        break if candidate.nil?

        current_log = log_likelihood(indicators, probabilities)
        previous_log = current_log
        if step.map { |value| value.abs * step_size }.max <= @tolerance
          converged = true
          break
        end
      end

      [probabilities, iterations, previous_log, converged]
    end

    def row_probability(row, probabilities)
      row.each_with_index.inject(0.0) do |total, (included, index)|
        total + (included ? probabilities[index] : 0.0)
      end
    end

    def solve_symmetric(matrix, vector)
      n = vector.size
      augmented = matrix.each_with_index.map { |row, index| row + [vector[index]] }

      n.times do |pivot|
        selected = (pivot...n).max_by { |row| augmented[row][pivot].abs }
        augmented[pivot], augmented[selected] = augmented[selected], augmented[pivot]
        pivot_value = augmented[pivot][pivot]
        raise ArgumentError, "singular active-set Hessian" if pivot_value.abs < 1e-18

        (pivot + 1...n).each do |row|
          factor = augmented[row][pivot] / pivot_value
          (pivot..n).each do |column|
            augmented[row][column] -= factor * augmented[pivot][column]
          end
        end
      end

      answer = Array.new(n, 0.0)
      (n - 1).downto(0) do |row|
        answer[row] = augmented[row][n]
        (row + 1...n).each { |column| answer[row] -= augmented[row][column] * answer[column] }
        answer[row] /= augmented[row][row]
      end
      answer
    end

    def kkt_residual(indicators, probabilities, active_indices = nil)
      scores = column_scores(indicators, probabilities)
      positive = active_indices || probabilities.each_index.select do |index|
        probabilities[index] > @zero_tolerance
      end
      return @kkt_tolerance if positive.empty?
      target = scores.values_at(*positive).sum / positive.size.to_f
      positive.map do |index|
        score = scores[index]
        if probabilities[index] > @zero_tolerance
          (score - target).abs
        else
          0.0
        end
      end.max
    end

    def column_scores(indicators, probabilities)
      scores = Array.new(probabilities.size, 0.0)
      indicators.each do |row|
        row_probability = 0.0
        row.each_with_index do |included, index|
          row_probability += probabilities[index] if included
        end
        row.each_with_index do |included, index|
          scores[index] += 1.0 / row_probability if included
        end
      end
      scores
    end

    def normalize_probabilities!(probabilities)
      total = probabilities.sum
      return probabilities if total.zero?
      probabilities.map! { |probability| probability / total }
    end
  end
end
