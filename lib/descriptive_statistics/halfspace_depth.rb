module DescriptiveStatistics

  # Two-dimensional halfspace (Tukey) depth and depth regions.
  #
  # The halfspace depth of a point +p+ relative to a sample is the minimum,
  # over all closed halfplanes containing +p+, of the fraction of the sample
  # inside the halfplane. Boundary convention: halfplanes are closed, so a
  # sample point lying exactly on the boundary line counts towards the
  # halfplane (this is the standard Tukey convention; consequently every
  # sample point has depth at least 1/n). Duplicate points count with their
  # multiplicity and collinear samples are supported.
  #
  # All computation is exact: coordinates are converted with +to_r+ and every
  # comparison is done in rational arithmetic, so integer or rational input
  # gives exact results.
  #
  #   points = [[0, 0], [4, 0], [4, 4], [0, 4]]
  #   DescriptiveStatistics::HalfspaceDepth.depth(points, [2, 2])  # => 0.5
  #   DescriptiveStatistics::HalfspaceDepth.depth_region(points, 2)
  #   # => #<DepthRegion type=:point vertices=[[2.0, 2.0]]>
  module HalfspaceDepth

    # The region of the plane where the halfspace depth is at least a given
    # threshold. +type+ is one of :empty, :point, :segment or :polygon;
    # degenerate regions are reported as :point or :segment rather than
    # being padded into zero-area polygons. +vertices+ holds the defining
    # points as [x, y] float pairs: none for :empty, one for :point, the two
    # endpoints for :segment, and the counter-clockwise corners for :polygon.
    class DepthRegion
      attr_reader :type, :vertices

      def initialize(type, vertices = [])
        @type = type
        @vertices = vertices
      end

      def empty?
        type == :empty
      end

      def point?
        type == :point
      end

      def segment?
        type == :segment
      end

      def polygon?
        type == :polygon
      end

      def ==(other)
        other.is_a?(DepthRegion) && type == other.type && vertices == other.vertices
      end

      alias eql? ==

      def hash
        [type, vertices].hash
      end

      def inspect
        "#<#{self.class.name} type=#{type.inspect} vertices=#{vertices.inspect}>"
      end
    end

    class << self

      # The Tukey halfspace depth of +point+ (an [x, y] pair) relative to
      # +points+ (an array of [x, y] pairs), as a Float in [0, 1]. Returns
      # nil for an empty sample.
      def depth(points, point)
        count = depth_count(points, point)
        return nil if count.nil?

        count.to_f / points.size
      end

      # The numerator of the halfspace depth: the minimum number of sample
      # points (counting duplicates) in any closed halfplane containing
      # +point+. Returns nil for an empty sample.
      def depth_count(points, point)
        sample = normalize_points(points)
        return nil if sample.empty?

        px, py = normalize_point(point)
        vectors = sample.map { |x, y| [x - px, y - py] }
        vectors.reject! { |vx, vy| vx.zero? && vy.zero? }
        coincident = sample.size - vectors.size
        return sample.size if vectors.empty?

        # For a vector v, the closed halfplane count is minimized just after
        # the boundary sweeps past v, i.e. by the open-left/closed-right
        # angular window (arg(v), arg(v) + pi].
        least = vectors.map do |v|
          vectors.count do |w|
            cross = v[0] * w[1] - v[1] * w[0]
            cross > 0 || (cross.zero? && v[0] * w[0] + v[1] * w[1] < 0)
          end
        end.min
        coincident + least
      end

      # The region of points whose halfspace depth is at least +k+. When +k+
      # is an Integer it is a depth count (1..n); when it is a Float or
      # Rational in (0, 1] it is a normalized depth, converted to a count by
      # multiplying by the sample size and rounding up. Returns a
      # DepthRegion; degenerate results come back as :point or :segment.
      def depth_region(points, k)
        sample = normalize_points(points)
        return DepthRegion.new(:empty) if sample.empty?

        count = threshold_count(k, sample.size)
        raise ArgumentError, "depth threshold must be positive" if count < 1
        return DepthRegion.new(:empty) if count > sample.size

        distinct = sample.uniq
        return DepthRegion.new(:point, floats([distinct.first])) if distinct.size == 1

        line = collinear_direction(sample)
        return collinear_region(sample, line, count) if line

        polygon_region(sample, distinct, count)
      end

      private

      def normalize_points(points)
        points.map { |point| normalize_point(point) }
      end

      def normalize_point(point)
        x, y = point
        [to_rational(x), to_rational(y)]
      end

      def to_rational(value)
        value.to_r
      rescue NoMethodError
        raise ArgumentError, "coordinates must be numeric, got #{value.inspect}"
      end

      def threshold_count(k, n)
        case k
        when Integer
          k
        when Float, Rational
          raise ArgumentError, "normalized depth must be in (0, 1]" if k <= 0 || k > 1

          (k * n).ceil
        else
          raise ArgumentError, "depth threshold must be an Integer count or a normalized depth in (0, 1]"
        end
      end

      def floats(vertices)
        vertices.map { |x, y| [x.to_f, y.to_f] }
      end

      # Returns a direction vector along the sample line when every point is
      # collinear, otherwise nil. Assumes at least two distinct points.
      def collinear_direction(sample)
        first = sample.first
        second = sample.find { |p| p != first }
        direction = [second[0] - first[0], second[1] - first[1]]
        collinear = sample.all? do |p|
          direction[0] * (p[1] - first[1]) == direction[1] * (p[0] - first[0])
        end
        collinear ? direction : nil
      end

      # Depth region of a collinear sample: a (possibly degenerate) segment
      # between the k-th and the (n - k + 1)-th order statistic along the line.
      def collinear_region(sample, direction, k)
        n = sample.size
        return DepthRegion.new(:empty) if k > n - k + 1

        origin = sample.first
        squared = direction[0] * direction[0] + direction[1] * direction[1]
        projections = sample.map do |p|
          (p[0] - origin[0]) * direction[0] + (p[1] - origin[1]) * direction[1]
        end.sort
        endpoints = [projections[k - 1], projections[n - k]].map do |t|
          [origin[0] + t * direction[0] / squared, origin[1] + t * direction[1] / squared]
        end
        if endpoints[0] == endpoints[1]
          DepthRegion.new(:point, floats([endpoints[0]]))
        else
          DepthRegion.new(:segment, floats(endpoints))
        end
      end

      # Depth region of a non-collinear sample: the intersection of every
      # closed halfplane that contains at least n - k + 1 sample points. It
      # suffices to consider halfplanes whose boundary line passes through
      # two distinct sample points.
      def polygon_region(sample, distinct, k)
        constraints = []
        distinct.combination(2) do |a, b|
          dx = b[0] - a[0]
          dy = b[1] - a[1]
          left = 0
          right = 0
          sample.each do |p|
            side = dx * (p[1] - a[1]) - dy * (p[0] - a[0])
            left += 1 if side > 0
            right += 1 if side < 0
          end
          constant = -dy * a[0] + dx * a[1]
          constraints << [-dy, dx, constant] if right <= k - 1
          constraints << [dy, -dx, -constant] if left <= k - 1
        end

        xs = sample.map(&:first)
        ys = sample.map(&:last)
        box = [[xs.min, ys.min], [xs.max, ys.min], [xs.max, ys.max], [xs.min, ys.max]]
        polygon = constraints.reduce(box) do |region, (nx, ny, constant)|
          clip(region, nx, ny, constant)
        end
        build_region(polygon)
      end

      # Keeps the part of a convex polygon inside the closed halfplane
      # nx*x + ny*y >= constant (Sutherland-Hodgman, exact arithmetic).
      def clip(polygon, nx, ny, constant)
        return [] if polygon.empty?

        clipped = []
        polygon.each_with_index do |start, index|
          finish = polygon[(index + 1) % polygon.size]
          ds = nx * start[0] + ny * start[1] - constant
          de = nx * finish[0] + ny * finish[1] - constant
          clipped << start if ds >= 0
          next if (ds >= 0) == (de >= 0)

          t = ds / (ds - de)
          clipped << [start[0] + t * (finish[0] - start[0]), start[1] + t * (finish[1] - start[1])]
        end
        clipped
      end

      # Classifies a clipped convex vertex list as empty, point, segment or
      # polygon, collapsing zero-area degeneracies instead of padding them.
      def build_region(vertices)
        pts = []
        vertices.each { |p| pts << p if pts.empty? || pts.last != p }
        pts.pop if pts.size > 1 && pts.first == pts.last

        while pts.size >= 3
          kept = []
          pts.each_with_index do |p, index|
            prev = pts[index - 1]
            nxt = pts[(index + 1) % pts.size]
            cross = (p[0] - prev[0]) * (nxt[1] - p[1]) - (p[1] - prev[1]) * (nxt[0] - p[0])
            between = (p[0] - prev[0]) * (p[0] - nxt[0]) + (p[1] - prev[1]) * (p[1] - nxt[1])
            kept << p unless cross.zero? && between <= 0
          end
          break if kept.size == pts.size

          pts = kept
        end

        case pts.size
        when 0 then return DepthRegion.new(:empty)
        when 1 then return DepthRegion.new(:point, floats(pts))
        when 2 then return DepthRegion.new(:segment, floats(pts))
        end

        area2 = pts.each_with_index.sum do |p, index|
          nxt = pts[(index + 1) % pts.size]
          p[0] * nxt[1] - p[1] * nxt[0]
        end
        if area2.zero?
          extremes = [pts.min, pts.max]
          return extremes[0] == extremes[1] ? DepthRegion.new(:point, floats([extremes[0]])) : DepthRegion.new(:segment, floats(extremes))
        end
        pts.reverse! if area2 < 0
        DepthRegion.new(:polygon, floats(pts))
      end

    end

  end

end
