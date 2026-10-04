module DescriptiveStatistics
  # Modified Band Depth (Lopez-Pintado & Romo, 2009) for a set of curves
  # observed on the same grid.
  #
  # Each curve is an array of values, one per grid position. Missing
  # observations are represented by nil and are never treated as zero:
  # a (grid position, curve pair) combination only enters a curve's
  # denominator when the curve itself and both curves of the pair are
  # observed at that position.
  #
  # The depth of a curve is the proportion of those valid
  # (position, pair) combinations in which the curve lies inside the
  # band (envelope) spanned by the two other curves of the pair.
  #
  #   curves = [[1, 2, 3], [3, 2, 1], [2, 2, nil]]
  #   mbd = DescriptiveStatistics::ModifiedBandDepth.new(curves)
  #   mbd.ranking            # => indices ordered from deepest outwards
  #   mbd.central_envelope(0.5) # => pointwise envelope of the deepest half
  class ModifiedBandDepth
    # Per-curve depth accounting. +inside+ and +total+ are the exact
    # numerator and denominator of the depth proportion, so it is always
    # explicit which (position, curve pair) combinations participated.
    class Depth
      attr_reader :index, :inside, :total

      def initialize(index, inside, total)
        @index = index
        @inside = inside
        @total = total
      end

      def value
        total.zero? ? nil : inside.to_f / total
      end
    end

    attr_reader :curves, :grid_size

    def initialize(curves)
      @curves = curves.map { |curve| curve.to_a.map { |value| value.nil? ? nil : Float(value) } }
      raise ArgumentError, "at least three curves are required" if @curves.size < 3

      @grid_size = @curves.first.size
      unless @curves.all? { |curve| curve.size == @grid_size }
        raise ArgumentError, "all curves must be observed on the same grid (equal length)"
      end
    end

    # Depth accounting for every curve, in the original curve order.
    def depths
      @depths ||= compute_depths
    end

    # Depth proportions in the original curve order (nil when a curve has
    # no valid (position, pair) combination).
    def depth_values
      depths.map(&:value)
    end

    # Curve indices ordered from deepest to most peripheral. Ties keep
    # the original curve order (stable); curves without any valid
    # combination sort last.
    def ranking
      depths
        .sort_by
        .with_index { |depth, order| [depth.value.nil? ? 1 : 0, -(depth.value || 0.0), order] }
        .map(&:index)
    end

    # Pointwise envelope of the deepest +proportion+ of the curves.
    # Returns a hash with the selected curve indices and the lower/upper
    # envelope per grid position. Positions where none of the selected
    # curves is observed are nil rather than zero.
    def central_envelope(proportion)
      unless proportion.is_a?(Numeric) && proportion > 0 && proportion <= 1
        raise ArgumentError, "proportion must be in (0, 1]"
      end

      selected = ranking.first((proportion * curves.size).ceil)
      lower = []
      upper = []
      grid_size.times do |position|
        observed = selected.map { |index| curves[index][position] }.compact
        lower << (observed.empty? ? nil : observed.min)
        upper << (observed.empty? ? nil : observed.max)
      end
      { indices: selected, lower: lower, upper: upper }
    end

    private

    def compute_depths
      curves.each_index.map do |index|
        inside = 0
        total = 0
        each_band_pair(index) do |lower_index, upper_index|
          grid_size.times do |position|
            value = curves[index][position]
            first = curves[lower_index][position]
            second = curves[upper_index][position]
            next if value.nil? || first.nil? || second.nil?

            total += 1
            low, high = first < second ? [first, second] : [second, first]
            inside += 1 if value >= low && value <= high
          end
        end
        Depth.new(index, inside, total)
      end
    end

    # Yields every unordered pair of curve indices not containing +excluded+.
    def each_band_pair(excluded)
      (0...curves.size).each do |first|
        next if first == excluded

        ((first + 1)...curves.size).each do |second|
          next if second == excluded

          yield first, second
        end
      end
    end
  end
end
