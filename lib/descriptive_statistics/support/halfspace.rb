require 'descriptive_statistics/depth_region.rb'

module DescriptiveStatistics
  module Support
    # Two-dimensional Tukey (halfspace) depth and depth regions.
    #
    # Boundary convention: depth is measured with closed halfplanes, so a
    # sample lying exactly on the boundary line of a halfplane is counted
    # as belonging to that halfplane. Duplicate and collinear samples are
    # supported.
    module Halfspace
      extend self

      RELATIVE_TOLERANCE = 1e-9
      ANGULAR_TOLERANCE = 1e-12

      def to_point(object)
        pair = object.to_a
        raise ArgumentError, "points must be two-element [x, y] pairs" unless pair.size == 2

        [Float(pair[0]), Float(pair[1])]
      rescue NoMethodError, TypeError
        raise ArgumentError, "points must be two-element [x, y] pairs"
      end

      def depth(points, query)
        coincident = 0
        angles = []
        points.each do |point|
          dx = point[0] - query[0]
          dy = point[1] - query[1]
          if dx.zero? && dy.zero?
            coincident += 1
          else
            angles << Math.atan2(dy, dx)
          end
        end
        return coincident if angles.empty?

        angles.sort!
        doubled = angles + angles.map { |angle| angle + 2.0 * Math::PI }
        smallest = angles.size
        angles.size.times do |index|
          first = doubled.bsearch_index { |angle| angle > doubled[index] + ANGULAR_TOLERANCE } || doubled.size
          last = doubled.bsearch_index { |angle| angle > doubled[index] + Math::PI + ANGULAR_TOLERANCE } || doubled.size
          count = last - first
          smallest = count if count < smallest
        end
        coincident + smallest
      end

      def region(points, k)
        return DepthRegion.new(:empty) if k > points.size

        distinct = points.uniq
        return DepthRegion.new(:point, [distinct.first]) if distinct.size == 1

        origin = distinct.first
        direction = subtract(distinct[1], origin)
        spread = tolerance(points) * magnitude(direction)
        if distinct.all? { |point| cross(direction, subtract(point, origin)).abs <= spread }
          return collinear_region(points, direction, k)
        end

        polygon = clip(bounding_box(points), halfplanes(points, k), tolerance(points))
        classify(polygon, points, k)
      end

      private

      def collinear_region(points, direction, k)
        return DepthRegion.new(:empty) if k - 1 > points.size - k

        sorted = points.sort_by { |point| dot(point, direction) }
        low = sorted[k - 1]
        high = sorted[points.size - k]
        return DepthRegion.new(:point, [low]) if low == high

        DepthRegion.new(:segment, [low, high])
      end

      def halfplanes(points, k)
        planes = []
        points.each do |a|
          points.each do |b|
            ux = a[1] - b[1]
            uy = b[0] - a[0]
            next if ux.zero? && uy.zero?

            projections = points.map { |point| dot(point, [ux, uy]) }.sort
            planes << [ux, uy, projections[k - 1]]
          end
        end
        planes
      end

      def clip(polygon, planes, tol)
        planes.each do |(ux, uy, c)|
          break if polygon.empty?

          margin = tol * magnitude([ux, uy])
          output = []
          polygon.each_with_index do |current, index|
            previous = polygon[index - 1]
            current_inside = dot(current, [ux, uy]) >= c - margin
            previous_inside = dot(previous, [ux, uy]) >= c - margin
            if current_inside
              output << intersection(previous, current, ux, uy, c) unless previous_inside
              output << current
            elsif previous_inside
              output << intersection(previous, current, ux, uy, c)
            end
          end
          polygon = output
        end
        polygon
      end

      def intersection(a, b, ux, uy, c)
        da = dot(a, [ux, uy]) - c
        db = dot(b, [ux, uy]) - c
        t = da / (da - db)
        [a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]
      end

      def classify(polygon, points, k)
        tol = tolerance(points)
        vertices = dedupe(polygon, tol)
        case vertices.size
        when 0
          DepthRegion.new(:empty)
        when 1
          point_region(vertices.first, points, k)
        when 2
          segment_region(vertices, points, k)
        else
          if signed_area(vertices).abs <= tol * scale(points)
            a, b = vertices.combination(2).max_by { |u, v| distance(u, v) }
            if distance(a, b) <= tol
              point_region(centroid(vertices), points, k)
            else
              segment_region([a, b], points, k)
            end
          else
            DepthRegion.new(:polygon, vertices)
          end
        end
      end

      def point_region(point, points, k)
        depth(points, point) >= k ? DepthRegion.new(:point, [point]) : DepthRegion.new(:empty)
      end

      def segment_region(endpoints, points, k)
        midpoint = centroid(endpoints)
        depth(points, midpoint) >= k ? DepthRegion.new(:segment, endpoints) : DepthRegion.new(:empty)
      end

      def dedupe(polygon, tol)
        vertices = polygon.each_with_object([]) do |vertex, kept|
          kept << vertex if kept.empty? || distance(vertex, kept.last) > tol
        end
        vertices.pop while vertices.size > 1 && distance(vertices.first, vertices.last) <= tol
        vertices
      end

      def bounding_box(points)
        xs = points.map { |point| point[0] }
        ys = points.map { |point| point[1] }
        margin = 10.0 * [xs.max - xs.min, ys.max - ys.min, 1.0].max
        minx = xs.min - margin
        maxx = xs.max + margin
        miny = ys.min - margin
        maxy = ys.max + margin
        [[minx, miny], [maxx, miny], [maxx, maxy], [minx, maxy]]
      end

      def centroid(vertices)
        [vertices.inject(0.0) { |acc, v| acc + v[0] } / vertices.size,
         vertices.inject(0.0) { |acc, v| acc + v[1] } / vertices.size]
      end

      def signed_area(vertices)
        total = vertices.each_with_index.inject(0.0) do |acc, (v, i)|
          acc + cross(v, vertices[(i + 1) % vertices.size])
        end
        total / 2.0
      end

      def scale(points)
        [1.0, points.map { |point| [point[0].abs, point[1].abs].max }.max].max
      end

      def tolerance(points)
        RELATIVE_TOLERANCE * scale(points)
      end

      def subtract(a, b)
        [a[0] - b[0], a[1] - b[1]]
      end

      def dot(a, b)
        a[0] * b[0] + a[1] * b[1]
      end

      def cross(a, b)
        a[0] * b[1] - a[1] * b[0]
      end

      def magnitude(a)
        Math.hypot(a[0], a[1])
      end

      def distance(a, b)
        magnitude(subtract(a, b))
      end
    end
  end
end
