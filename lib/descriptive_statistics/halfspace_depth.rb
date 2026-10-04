module DescriptiveStatistics
  # Tukey (halfspace) depth of +point+ within a collection of [x, y] pairs.
  #
  # The depth is the minimum number of samples contained in any closed
  # halfplane that contains +point+: samples lying exactly on the boundary
  # line of a halfplane are counted as inside it. Duplicate and collinear
  # samples are supported. Returns an Integer; points outside the convex
  # hull of the samples have depth 0.
  def halfspace_depth(point, collection = self, &block)
    values = Support::extract(collection, &block).to_a
    return DescriptiveStatistics.halfspace_depth_empty_collection_default_value if values.empty?

    query = Support::Halfspace.to_point(point)
    points = values.map { |value| Support::Halfspace.to_point(value) }
    Support::Halfspace.depth(points, query)
  end
end
