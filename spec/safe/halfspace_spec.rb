require 'safe/spec_helper'

describe "Halfspace depth and depth regions" do
  require 'descriptive_statistics/safe'

  # Independent reference: the Tukey depth is the minimum, over all closed
  # halfplanes containing the query, of the number of samples inside. The
  # minimum is attained at a critical direction (normal to some sample
  # vector) or immediately past one, so enumerating those directions
  # suffices. Integer arithmetic keeps the enumeration exact.
  def brute_force_depth(points, query)
    vectors = points.map { |p| [p[0] - query[0], p[1] - query[1]] }
    nonzero = vectors.reject { |v| v[0] == 0 && v[1] == 0 }
    return points.size if nonzero.empty?

    normals = nonzero.flat_map { |v| [[-v[1], v[0]], [v[1], -v[0]]] }
    normals.flat_map do |u|
      closed = vectors.count { |w| w[0] * u[0] + w[1] * u[1] >= 0 }
      past = vectors.count do |w|
        dot = w[0] * u[0] + w[1] * u[1]
        dot > 0 || (dot == 0 && u[0] * w[1] - u[1] * w[0] >= 0)
      end
      [closed, past]
    end.min
  end

  def random_affine
    begin
      a, b, c, d = Array.new(4) { rand(-2.0..2.0) }
    end until (a * d - b * c).abs > 0.5
    t = [rand(-3.0..3.0), rand(-3.0..3.0)]
    ->(p) { [a * p[0] + b * p[1] + t[0], c * p[0] + d * p[1] + t[1]] }
  end

  let(:square) { [[0, 0], [1, 0], [1, 1], [0, 1]] }

  context "halfspace_depth" do
    it "computes known values on a unit square" do
      data = DescriptiveStatistics::Stats.new(square)
      expect(data.halfspace_depth([0.5, 0.5])).to eql(2)
      expect(data.halfspace_depth([0, 0])).to eql(1)
      expect(data.halfspace_depth([2, 2])).to eql(0)
    end

    it "counts samples on the boundary line as inside the halfplane" do
      data = DescriptiveStatistics::Stats.new([[0, 0], [1, 0], [2, 0]])
      expect(data.halfspace_depth([1, 0])).to eql(2)
      expect(data.halfspace_depth([0, 0])).to eql(1)
      expect(data.halfspace_depth([0.5, 0])).to eql(1)
    end

    it "supports duplicate samples" do
      data = DescriptiveStatistics::Stats.new([[0, 0], [0, 0], [2, 0]])
      expect(data.halfspace_depth([0, 0])).to eql(2)
      expect(data.halfspace_depth([1, 0])).to eql(1)
    end

    it "is available as a module function" do
      expect(DescriptiveStatistics.halfspace_depth([0.5, 0.5], square)).to eql(2)
    end

    it "maps objects with a block" do
      records = [{ :p => [0, 0] }, { :p => [1, 0] }, { :p => [1, 1] }, { :p => [0, 1] }]
      expect(DescriptiveStatistics::Stats.new(records).halfspace_depth([0.5, 0.5]) { |r| r[:p] }).to eql(2)
    end

    it "returns the empty collection default value" do
      expect(DescriptiveStatistics::Stats.new([]).halfspace_depth([0, 0])).to be_nil
    end

    it "rejects malformed points" do
      expect { DescriptiveStatistics::Stats.new([[1, 2], [3]]).halfspace_depth([0, 0]) }.to raise_error(ArgumentError)
      expect { DescriptiveStatistics::Stats.new([[1, 2]]).halfspace_depth(5) }.to raise_error(ArgumentError)
    end

    it "matches brute-force enumeration of critical directions on small samples" do
      srand(7)
      200.times do
        points = Array.new(rand(2..10)) { [rand(-4..4), rand(-4..4)] }
        points << points.sample if rand < 0.4
        if rand < 0.3
          rand(2..4).times { |x| points << [x, 2 * x + 1] }
        end
        query = [rand(-5..5), rand(-5..5)]
        expect(DescriptiveStatistics::Stats.new(points).halfspace_depth(query))
          .to eql(brute_force_depth(points, query))
      end
    end

    it "is invariant under affine transformations" do
      srand(99)
      100.times do
        points = Array.new(rand(3..10)) { [rand(-5.0..5.0), rand(-5.0..5.0)] }
        query = [rand(-5.0..5.0), rand(-5.0..5.0)]
        map = random_affine
        expect(DescriptiveStatistics::Stats.new(points.map(&map)).halfspace_depth(map.call(query)))
          .to eql(DescriptiveStatistics::Stats.new(points).halfspace_depth(query))
      end
    end
  end

  context "depth_region" do
    it "returns the convex hull as a polygon for k = 1" do
      region = DescriptiveStatistics::Stats.new(square).depth_region(1)
      expect(region.kind).to eql(:polygon)
      expect(region.polygon?).to be(true)
      expect(region.vertices.size).to eql(4)
      square.each do |corner|
        expect(region.vertices.any? { |v| (v[0] - corner[0]).abs < 1e-9 && (v[1] - corner[1]).abs < 1e-9 }).to be(true)
      end
    end

    it "returns a point when the deepest region is a single location" do
      region = DescriptiveStatistics::Stats.new(square).depth_region(2)
      expect(region.kind).to eql(:point)
      expect(region.point?).to be(true)
      expect(region.vertices.size).to eql(1)
      expect(region.vertices.first[0]).to be_within(1e-9).of(0.5)
      expect(region.vertices.first[1]).to be_within(1e-9).of(0.5)
    end

    it "returns a segment for collinear samples instead of a flat polygon" do
      region = DescriptiveStatistics::Stats.new([[0, 0], [1, 0], [2, 0], [3, 0]]).depth_region(2)
      expect(region.kind).to eql(:segment)
      expect(region.segment?).to be(true)
      expect(region.vertices).to eql([[1.0, 0.0], [2.0, 0.0]])
    end

    it "returns a point for the median of an odd collinear sample" do
      region = DescriptiveStatistics::Stats.new([[0, 0], [1, 0], [2, 0], [3, 0], [4, 0]]).depth_region(3)
      expect(region.kind).to eql(:point)
      expect(region.vertices).to eql([[2.0, 0.0]])
    end

    it "returns a point when all samples coincide" do
      region = DescriptiveStatistics::Stats.new([[1, 1], [1, 1], [1, 1]]).depth_region(2)
      expect(region.kind).to eql(:point)
      expect(region.vertices).to eql([[1.0, 1.0]])
    end

    it "returns an empty region when k exceeds the maximum depth" do
      expect(DescriptiveStatistics::Stats.new(square).depth_region(3).kind).to eql(:empty)
      expect(DescriptiveStatistics::Stats.new(square).depth_region(3).empty?).to be(true)
      expect(DescriptiveStatistics::Stats.new([[0, 0], [1, 0], [2, 0], [3, 0]]).depth_region(3).kind).to eql(:empty)
    end

    it "rejects non-positive or non-integer levels" do
      data = DescriptiveStatistics::Stats.new(square)
      expect { data.depth_region(0) }.to raise_error(ArgumentError)
      expect { data.depth_region(1.5) }.to raise_error(ArgumentError)
    end

    it "returns the empty collection default value" do
      expect(DescriptiveStatistics::Stats.new([]).depth_region(1)).to be_nil
    end

    it "has depth at least k inside the region and less just outside" do
      points = [[0, 0], [4, 0], [0, 4], [4, 4], [1, 1], [1, 3], [3, 1], [3, 3]]
      data = DescriptiveStatistics::Stats.new(points)
      region = data.depth_region(2)
      expect(region.kind).to eql(:polygon)
      cx = region.vertices.sum { |v| v[0] } / region.vertices.size
      cy = region.vertices.sum { |v| v[1] } / region.vertices.size
      region.vertices.each do |v|
        inner = [cx + 0.99 * (v[0] - cx), cy + 0.99 * (v[1] - cy)]
        outer = [cx + 1.5 * (v[0] - cx), cy + 1.5 * (v[1] - cy)]
        expect(data.halfspace_depth(inner)).to be >= 2
        expect(data.halfspace_depth(outer)).to be < 2
      end
    end

    it "is equivariant under affine transformations" do
      srand(5)
      50.times do
        points = Array.new(9) { [rand(-5.0..5.0), rand(-5.0..5.0)] }
        map = random_affine
        k = [1, 2, 3].sample
        region = DescriptiveStatistics::Stats.new(points).depth_region(k)
        mapped_region = DescriptiveStatistics::Stats.new(points.map(&map)).depth_region(k)
        expect(mapped_region.kind).to eql(region.kind)
        expect(mapped_region.vertices.size).to eql(region.vertices.size)
        region.vertices.map(&map).each do |m|
          expect(mapped_region.vertices.any? { |v| (v[0] - m[0]).abs < 1e-6 && (v[1] - m[1]).abs < 1e-6 }).to be(true)
        end
      end
    end
  end
end
