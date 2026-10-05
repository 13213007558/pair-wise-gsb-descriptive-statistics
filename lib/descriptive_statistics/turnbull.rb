module DescriptiveStatistics

  # Turnbull's nonparametric maximum likelihood estimator (NPMLE) for the
  # distribution of a random variable whose observations may be exact,
  # left censored, right censored or interval censored.
  #
  # Every observation is an interval with explicitly declared open/closed
  # endpoints:
  #
  #   observations = [
  #     DescriptiveStatistics::Turnbull.exact(2.0),
  #     DescriptiveStatistics::Turnbull.interval(1.0, 3.0, left_closed: false, right_closed: true),
  #     DescriptiveStatistics::Turnbull.left_censored(1.5),
  #     DescriptiveStatistics::Turnbull.right_censored(4.0, closed: false)
  #   ]
  #
  #   fit = DescriptiveStatistics::Turnbull.new(observations).fit
  #   fit.support       # => [#<struct interval=(1.0, 2.0], mass=0.5>, ...]
  #   fit.log_likelihood
  #   fit.diagnostics   # => { iterations: ..., converged: ..., ... }
  #   fit.cdf(2.5)      # => 0.25..0.75
  #
  # The estimated probability mass is reported per maximal support interval.
  # The likelihood only identifies the total mass of each such interval, not
  # its placement inside, so the mass is never collapsed onto an arbitrary
  # point (such as the midpoint). Consequently +cdf+ returns a range of
  # achievable values whenever the query point falls inside a support
  # interval.
  class Turnbull

    # An interval of the real line with explicitly declared open/closed
    # endpoints. Infinite endpoints are always open.
    class Interval
      attr_reader :left, :right

      def initialize(left, right, left_closed: true, right_closed: true)
        left = Float(left)
        right = Float(right)
        raise ArgumentError, "endpoints must not be NaN" if left.nan? || right.nan?
        raise ArgumentError, "left endpoint must not exceed right endpoint" if left > right
        if left == right && !(left_closed && right_closed)
          raise ArgumentError, "a degenerate interval must be closed on both sides"
        end
        if left == -Float::INFINITY && left_closed
          raise ArgumentError, "an infinite left endpoint cannot be closed"
        end
        if right == Float::INFINITY && right_closed
          raise ArgumentError, "an infinite right endpoint cannot be closed"
        end
        @left = left
        @right = right
        @left_closed = left_closed
        @right_closed = right_closed
      end

      def left_closed?
        @left_closed
      end

      def right_closed?
        @right_closed
      end

      def ==(other)
        other.is_a?(Interval) &&
          left == other.left && right == other.right &&
          left_closed? == other.left_closed? && right_closed? == other.right_closed?
      end
      alias eql? ==

      def hash
        [left, right, left_closed?, right_closed?].hash
      end

      def to_s
        "#{left_closed? ? '[' : '('}#{left}, #{right}#{right_closed? ? ']' : ')'}"
      end
      alias inspect to_s
    end

    # A maximal support interval of the NPMLE together with its total
    # probability mass. The placement of the mass inside the interval is not
    # identifiable from the data.
    SupportInterval = Struct.new(:interval, :mass)

    # Probability masses below this threshold are dropped from the support.
    MASS_THRESHOLD = 1e-8

    class << self
      # A general interval observation, e.g. interval(1.0, 3.0,
      # left_closed: false, right_closed: true) for (1.0, 3.0].
      def interval(left, right, left_closed: true, right_closed: true)
        Interval.new(left, right, left_closed: left_closed, right_closed: right_closed)
      end

      # An exact (uncensored) observation of +value+.
      def exact(value)
        Interval.new(value, value, left_closed: true, right_closed: true)
      end

      # A left censored observation: the event occurred at or before
      # (+closed: true+) or strictly before (+closed: false+) +upper+.
      def left_censored(upper, closed: true)
        Interval.new(-Float::INFINITY, upper, left_closed: false, right_closed: closed)
      end

      # A right censored observation: the event occurred at or after
      # (+closed: true+) or strictly after (+closed: false+) +lower+.
      def right_censored(lower, closed: true)
        Interval.new(lower, Float::INFINITY, left_closed: closed, right_closed: false)
      end
    end

    attr_reader :observations, :tolerance, :max_iterations

    # +observations+ is a collection of Turnbull::Interval objects (usually
    # built via the Turnbull.interval/exact/left_censored/right_censored
    # constructors). Plain numerics are treated as exact observations.
    def initialize(observations, tolerance: 1e-9, max_iterations: 10_000)
      observations = observations.to_a
      raise ArgumentError, "at least one observation is required" if observations.empty?
      @observations = observations.map { |obs| obs.is_a?(Interval) ? obs : Turnbull.exact(obs) }
      @tolerance = tolerance
      @max_iterations = max_iterations
      @fitted = false
    end

    # Runs the EM (self-consistency) iterations and returns self.
    def fit
      sets = elementary_sets
      membership = @observations.map { |obs| covered_indices(obs, sets) }
      n = @observations.size
      k = sets.size
      probabilities = Array.new(k, 1.0 / k)
      trace = []
      converged = false
      iterations = 0
      max_change = nil

      while iterations < @max_iterations && !converged
        updated = Array.new(k, 0.0)
        log_likelihood = 0.0
        membership.each do |indices|
          denom = indices.sum { |i| probabilities[i] }
          log_likelihood += Math.log(denom)
          indices.each { |i| updated[i] += probabilities[i] / denom }
        end
        updated.map! { |u| u / n }
        max_change = updated.zip(probabilities).map { |a, b| (a - b).abs }.max
        probabilities = updated
        trace << log_likelihood
        iterations += 1
        converged = max_change < @tolerance
      end

      @log_likelihood = trace.last
      @support = build_support(sets, probabilities, membership)
      @diagnostics = {
        iterations: iterations,
        converged: converged,
        log_likelihood: @log_likelihood,
        log_likelihood_trace: trace,
        max_probability_change: max_change,
        tolerance: @tolerance,
        max_iterations: @max_iterations
      }
      @fitted = true
      self
    end

    # The maximal support intervals and their total probability masses.
    def support
      ensure_fitted
      @support
    end

    def log_likelihood
      ensure_fitted
      @log_likelihood
    end

    def diagnostics
      ensure_fitted
      @diagnostics
    end

    # The range of values the estimated CDF can take at +x+. Both bounds
    # coincide unless +x+ falls strictly inside a support interval, in which
    # case the placement of that interval's mass is not identifiable.
    def cdf(x)
      cdf_bounds(x).first..cdf_bounds(x).last
    end

    def cdf_lower(x)
      cdf_bounds(x).first
    end

    def cdf_upper(x)
      cdf_bounds(x).last
    end

    private

    def ensure_fitted
      fit unless @fitted
    end

    # Partitions the real line into the elementary sets induced by the
    # finite observation endpoints: the open ray below the smallest
    # endpoint, each endpoint as a singleton, each open gap between
    # consecutive endpoints, and the open ray above the largest endpoint.
    def elementary_sets
      endpoints = @observations
        .flat_map { |obs| [obs.left, obs.right] }
        .reject { |e| e.infinite? }
        .uniq
        .sort
      return [Interval.new(-Float::INFINITY, Float::INFINITY, left_closed: false, right_closed: false)] if endpoints.empty?

      sets = [Interval.new(-Float::INFINITY, endpoints.first, left_closed: false, right_closed: false)]
      endpoints.each_with_index do |t, j|
        sets << Interval.new(t, t, left_closed: true, right_closed: true)
        if j + 1 < endpoints.size
          sets << Interval.new(t, endpoints[j + 1], left_closed: false, right_closed: false)
        end
      end
      sets << Interval.new(endpoints.last, Float::INFINITY, left_closed: false, right_closed: false)
      sets
    end

    def covered_indices(obs, sets)
      sets.each_index.select { |i| covers?(obs, sets[i]) }
    end

    def covers?(obs, set)
      if set.left == set.right
        t = set.left
        (obs.left < t || (obs.left == t && obs.left_closed?)) &&
          (obs.right > t || (obs.right == t && obs.right_closed?))
      else
        obs.left <= set.left && obs.right >= set.right
      end
    end

    # Groups elementary sets into maximal support intervals. Two adjacent
    # elementary sets belong to the same support interval exactly when the
    # same observations cover them, because then the likelihood cannot
    # distinguish mass placed on one from mass placed on the other. Only the
    # total mass of each group is identifiable. Masses are renormalized to
    # sum to one after negligible masses are dropped.
    def build_support(sets, probabilities, membership)
      patterns = Array.new(sets.size) { [] }
      membership.each_with_index do |indices, obs_index|
        indices.each { |i| patterns[i] << obs_index }
      end

      groups = []
      sets.each_with_index do |set, i|
        last = groups.last
        if last && last[:pattern] == patterns[i]
          last[:mass] += probabilities[i]
          last[:right] = set.right
          last[:right_closed] = set.right_closed?
        else
          groups << {
            pattern: patterns[i],
            mass: probabilities[i],
            left: set.left,
            left_closed: set.left_closed?,
            right: set.right,
            right_closed: set.right_closed?
          }
        end
      end

      groups.select! { |group| group[:mass] > MASS_THRESHOLD }
      total = groups.sum { |group| group[:mass] }
      groups.map do |group|
        interval = Interval.new(group[:left], group[:right],
          left_closed: group[:left_closed], right_closed: group[:right_closed])
        SupportInterval.new(interval, group[:mass] / total)
      end
    end

    def cdf_bounds(x)
      ensure_fitted
      lower = 0.0
      upper = 0.0
      @support.each do |entry|
        interval = entry.interval
        if interval.right <= x
          lower += entry.mass
          upper += entry.mass
        elsif interval.left > x || (interval.left == x && !interval.left_closed?)
          # no point of the interval is at or below x
        else
          upper += entry.mass
        end
      end
      [lower, upper]
    end

  end

end
