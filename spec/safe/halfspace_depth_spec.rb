require 'safe/spec_helper'

describe "DescriptiveStatistics::HalfspaceDepth" do
  require 'descriptive_statistics/safe'

  let(:depth) { DescriptiveStatistics::HalfspaceDepth }

  # Independent reference implementation used to cross-check the library:
  # it enumerates every critical direction (normals perpendicular to a
  # sample vector) and evaluates the closed-halfplane sample count just
  # before and just after each of them, using exact rational arithmetic.
  def brute_force_depth_count(points, query)
    vectors = points.map { |x, y| [x - query[0], y - query[1]] }
    coincident = vectors.count { |vx, vy| vx.zero? && vy.zero? }
    rest = vectors.reject { |vx, vy| vx.zero? && vy.zero? }
    return points.size if rest.empty?

    best = rest.size
    rest.each do |vx, vy|
      [[-vy, vx], [vy, -vx]].each do |ux, uy|
        [1, -1].each do |nudge|
          count = rest.count do |wx, wy|
            dot = ux * wx + uy * wy
            dot > 0 || (dot.zero? && nudge * (vx * wx + vy * wy) >= 0)
          end
          best = count if count < best
        end
      end
    end
    coincident + best
  end

  def region_contains?(region, x, y, eps = 1e-6)
    case region.type
    when :empty
      false
    when :point
      px, py = region.vertices.first
      (px - x).abs <= eps && (py - y).abs <= eps
    when :segment
      (ax, ay), (bx, by) = region.vertices
      cross = (bx - ax) * (y - ay) - (by - ay) * (x - ax)
      scale = [(bx - ax).abs + (by - ay).abs, 1.0].max
      cross.abs <= eps * scale &&
        x >= [ax, bx].min - eps && x <= [ax, bx].max + eps &&
        y >= [ay, by].min - eps && y <= [ay, by].max + eps
    when :polygon
      vertices = region.vertices
      vertices.each_with_index.all? do |(ax, ay), index|
        bx, by = vertices[(index + 1) % vertices.size]
        (bx - ax) * (y - ay) - (by - ay) * (x - ax) >= -eps
      end
    end
  end

  context "boundary convention" do
    it "counts samples on the boundary line towards the closed halfplane" do
      points = [[0, 0], [1, 0], [2, 0]]
      expect(depth.depth_count(points, [1, 0])).to eql(2)
      expect(depth.depth(points, [1, 0])).to eql(2.0 / 3)
    end

    it "gives every sample point depth at least 1/n" do
      points = [[0, 0], [4, 0], [4, 4], [0, 4]]
      points.each do |point|
        expect(depth.depth_count(points, point)).to be >= 1
      end
      expect(depth.depth(points, [0, 0])).to eql(0.25)
    end

    it "gives depth zero to points separable from the whole sample" do
      points = [[0, 0], [4, 0], [4, 4], [0, 4]]
      expect(depth.depth(points, [10, 10])).to eql(0.0)
    end

    it "returns nil for an empty sample" do
      expect(depth.depth([], [0, 0])).to be_nil
      expect(depth.depth_count([], [0, 0])).to be_nil
    end
  end

  context "with duplicate and collinear samples" do
    it "counts duplicates with their multiplicity" do
      expect(depth.depth([[0, 0], [0, 0], [1, 0]], [0, 0])).to eql(2.0 / 3)
      expect(depth.depth([[0, 0], [0, 0], [1, 0]], [1, 0])).to eql(1.0 / 3)
    end

    it "handles a fully coincident sample" do
      points = [[5, 5], [5, 5], [5, 5]]
      expect(depth.depth(points, [5, 5])).to eql(1.0)
      expect(depth.depth(points, [0, 0])).to eql(0.0)
    end

    it "handles collinear samples" do
      points = [[0, 0], [1, 0], [2, 0], [3, 0]]
      expect(depth.depth(points, [1, 0])).to eql(0.5)
      expect(depth.depth(points, [2, 0])).to eql(0.5)
      expect(depth.depth(points, [10, 0])).to eql(0.0)
    end
  end

  context "cross-checked against enumeration of critical directions" do
    it "matches the brute force count on small random samples" do
      rng = Random.new(20261004)
      200.times do
        n = rng.rand(1..8)
        points = Array.new(n) { [rng.rand(-3..3), rng.rand(-3..3)] }
        query = [rng.rand(-4..4), rng.rand(-4..4)]
        expect(depth.depth_count(points, query)).to eql(brute_force_depth_count(points, query))
      end
    end

    it "matches the brute force count on collinear and duplicate-heavy samples" do
      rng = Random.new(42)
      100.times do
        pool = [[0, 0], [1, 1], [2, 2], [3, 3], [0, 0], [1, 1]]
        points = Array.new(rng.rand(1..8)) { pool[rng.rand(pool.size)] }
        query = [rng.rand(-1..4), rng.rand(-1..4)]
        expect(depth.depth_count(points, query)).to eql(brute_force_depth_count(points, query))
      end
    end
  end

  context "affine invariance" do
    transforms = [
      [[2, 1], [1, 1], [3, -2]],
      [[1, 2], [0, 1], [0, 0]],
      [[3, 1], [5, 2], [-1, 4]],
      [[-2, 0], [0, 3], [5, 5]],
      [[0, -1], [1, 0], [2, 2]]
    ]

    def apply_transform(transform, point)
      (a, b), (c, d), (tx, ty) = transform
      x, y = point
      [a * x + b * y + tx, c * x + d * y + ty]
    end

    it "keeps the depth count unchanged under nonsingular affine maps" do
      rng = Random.new(7)
      samples = [
        [[0, 0], [4, 0], [4, 4], [0, 4]],
        [[0, 0], [1, 0], [2, 0], [3, 0]],
        [[0, 0], [0, 0], [1, 2], [2, 1], [3, 3]],
        Array.new(9) { [rng.rand(-3..3), rng.rand(-3..3)] }
      ]
      samples.each do |points|
        queries = points.first(2) + [[1, 1], [rng.rand(-2..2), rng.rand(-2..2)]]
        transforms.each do |transform|
          mapped = points.map { |point| apply_transform(transform, point) }
          queries.each do |query|
            expect(depth.depth_count(mapped, apply_transform(transform, query)))
              .to eql(depth.depth_count(points, query))
          end
        end
      end
    end
  end

  context "depth regions" do
    let(:square) { [[0, 0], [4, 0], [4, 4], [0, 4]] }

    it "returns the convex hull as a polygon for k = 1" do
      region = depth.depth_region(square, 1)
      expect(region.type).to eql(:polygon)
      expect(region.vertices).to eql([[0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0]])
    end

    it "returns a point, not a zero-area polygon, when the region is a singleton" do
      region = depth.depth_region(square, 2)
      expect(region.type).to eql(:point)
      expect(region.vertices).to eql([[2.0, 2.0]])
    end

    it "returns a segment, not a zero-area polygon, for collinear samples" do
      region = depth.depth_region([[0, 0], [1, 0], [2, 0], [3, 0]], 2)
      expect(region.type).to eql(:segment)
      expect(region.vertices).to eql([[1.0, 0.0], [2.0, 0.0]])
    end

    it "returns a segment hull for collinear samples at k = 1" do
      region = depth.depth_region([[0, 0], [3, 0], [1, 0], [2, 0]], 1)
      expect(region.type).to eql(:segment)
      expect(region.vertices).to eql([[0.0, 0.0], [3.0, 0.0]])
    end

    it "collapses to a point at the median of a collinear sample" do
      region = depth.depth_region([[0, 0], [1, 0], [2, 0], [3, 0], [4, 0]], 3)
      expect(region.type).to eql(:point)
      expect(region.vertices).to eql([[2.0, 0.0]])
    end

    it "returns a point for a fully coincident sample" do
      region = depth.depth_region([[5, 5], [5, 5], [5, 5]], 2)
      expect(region.type).to eql(:point)
      expect(region.vertices).to eql([[5.0, 5.0]])
    end

    it "returns an empty region when k exceeds the maximum depth" do
      expect(depth.depth_region(square, 3).type).to eql(:empty)
      expect(depth.depth_region([[0, 0], [6, 0], [0, 6]], 2).type).to eql(:empty)
      expect(depth.depth_region([], 1).type).to eql(:empty)
    end

    it "accepts a normalized depth threshold" do
      expect(depth.depth_region(square, 0.5)).to eql(depth.depth_region(square, 2))
      expect(depth.depth_region(square, 0.25)).to eql(depth.depth_region(square, 1))
    end

    it "rejects invalid thresholds" do
      expect { depth.depth_region(square, 0) }.to raise_error(ArgumentError)
      expect { depth.depth_region(square, -1) }.to raise_error(ArgumentError)
      expect { depth.depth_region(square, 1.5) }.to raise_error(ArgumentError)
      expect { depth.depth_region(square, "2") }.to raise_error(ArgumentError)
    end

    it "keeps region vertices at depth at least k" do
      [1, 2].each do |k|
        region = depth.depth_region(square, k)
        region.vertices.each do |vertex|
          expect(depth.depth_count(square, vertex)).to be >= k
        end
      end
    end

    it "matches depth-based membership on a grid of probe points" do
      rng = Random.new(99)
      samples = [
        [[0, 0], [4, 0], [4, 4], [0, 4]],
        [[0, 0], [6, 0], [0, 6]],
        [[0, 0], [4, 0], [4, 4], [0, 4], [2, 2]],
        [[0, 0], [1, 0], [2, 0], [3, 0]],
        [[1, 1], [1, 1], [3, 1], [1, 3], [3, 3]],
        Array.new(7) { [rng.rand(0..5), rng.rand(0..5)] },
        Array.new(9) { [rng.rand(0..6), rng.rand(0..6)] }
      ]
      samples.each do |points|
        (1..3).each do |k|
          region = depth.depth_region(points, k)
          (-2..14).each do |i|
            (-2..14).each do |j|
              x = i / 2.0
              y = j / 2.0
              expect(region_contains?(region, x, y)).to eql(depth.depth_count(points, [x, y]) >= k)
            end
          end
        end
      end
    end
  end
end
