require 'monkeypatch/spec_helper'

describe "DescriptiveStatistics" do
  require 'descriptive_statistics'

  context "with an array of two-dimensional points" do
    subject { [[0, 0], [1, 0], [1, 1], [0, 1]] }

    it "computes the halfspace depth of a query point" do
      expect(subject.halfspace_depth([0.5, 0.5])).to eql(2)
      expect(subject.halfspace_depth([2, 2])).to eql(0)
    end

    it "computes depth regions" do
      expect(subject.depth_region(1).kind).to eql(:polygon)
      expect(subject.depth_region(2).kind).to eql(:point)
      expect(subject.depth_region(3).kind).to eql(:empty)
    end

    it "does not force degenerate regions into polygons" do
      expect([[0, 0], [1, 0], [2, 0], [3, 0]].depth_region(2).kind).to eql(:segment)
      expect([[0, 0], [1, 0], [2, 0], [3, 0], [4, 0]].depth_region(3).kind).to eql(:point)
    end
  end
end
