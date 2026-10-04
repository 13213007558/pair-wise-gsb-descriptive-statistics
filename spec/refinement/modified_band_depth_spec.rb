require 'refinement/spec_helper'

describe "DescriptiveStatistics modified band depth refinements" do
  require 'descriptive_statistics/refinement'

  class ModifiedBandDepthCalculator
    using DescriptiveStatistics::Refinement.new(Array)

    def self.depth(curves)
      curves.modified_band_depth
    end

    def self.order(curves)
      curves.modified_band_depth_order
    end

    def self.envelope(curves)
      curves.central_envelope
    end
  end

  it "calculates depth through a refinement" do
    curves = [[-1, 4, -1], [4, -1, 4], [1, 1, 1]]

    expect(ModifiedBandDepthCalculator.depth(curves)).to eq([0.0, 0.0, 1.0])
    expect(ModifiedBandDepthCalculator.order(curves)).to eq([2, 0, 1])
  end

  it "calculates the central envelope through a refinement" do
    curves = [[-1, 4, -1], [4, -1, 4], [1, 1, 1]]

    expect(ModifiedBandDepthCalculator.envelope(curves)).to eq([
      [-1, 1, -1], [1, 4, 1]
    ])
  end
end
