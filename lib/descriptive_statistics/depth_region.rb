module DescriptiveStatistics
  # The set of points whose halfspace (Tukey) depth is at least k.
  #
  # The region is convex but may degenerate: +kind+ is one of :polygon,
  # :segment, :point or :empty, and +vertices+ holds the boundary vertices
  # (counter-clockwise for polygons, the two endpoints for segments, the
  # single location for points). Degenerate regions are reported as
  # segments or points rather than zero-area polygons.
  class DepthRegion
    attr_reader :kind, :vertices

    def initialize(kind, vertices = [])
      @kind = kind
      @vertices = vertices
    end

    def empty?
      kind == :empty
    end

    def point?
      kind == :point
    end

    def segment?
      kind == :segment
    end

    def polygon?
      kind == :polygon
    end

    def to_a
      vertices
    end

    def ==(other)
      other.is_a?(DepthRegion) && kind == other.kind && vertices == other.vertices
    end
    alias_method :eql?, :==

    def hash
      [kind, vertices].hash
    end

    def inspect
      "#<DescriptiveStatistics::DepthRegion #{kind} #{vertices.inspect}>"
    end
    alias_method :to_s, :inspect
  end

  def depth_region(k, collection = self, &block)
    values = Support::extract(collection, &block).to_a
    return DescriptiveStatistics.depth_region_empty_collection_default_value if values.empty?

    level = k.is_a?(Float) && k == k.to_i ? k.to_i : k
    raise ArgumentError, "depth level must be a positive integer" unless level.is_a?(Integer) && level >= 1

    points = values.map { |value| Support::Halfspace.to_point(value) }
    Support::Halfspace.region(points, level)
  end
end
