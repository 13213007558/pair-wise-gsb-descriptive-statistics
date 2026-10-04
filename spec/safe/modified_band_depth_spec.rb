require 'safe/spec_helper'

describe "DescriptiveStatistics modified band depth" do
  require 'descriptive_statistics/safe'

  let(:curves) do
    [
      [-1, 4, -1],
      [4, -1, 4],
      [1, 1, 1]
    ]
  end

  it "calculates depth through module functions" do
    expect(DescriptiveStatistics.modified_band_depth(curves)).to eq(
      [0.0, 0.0, 1.0]
    )
  end

  it "returns the central envelope through module functions" do
    expect(DescriptiveStatistics.central_envelope(curves, 0.5)).to eq([
      [-1, 1, -1], [1, 4, 1]
    ])
  end

  it "exposes denominator details through module functions" do
    curves = [[1, nil, 1], [1, 0, 1], [1, 2, 1]]
    details = DescriptiveStatistics.modified_band_depth_details(curves)

    expect(details.first.fetch(:depth)).to eq(1.0)
    expect(details.first.fetch(:bands).first.fetch(:observed_positions)).to eq([0, 2])
    expect(details.first.fetch(:bands).first.fetch(:contained_positions)).to eq([0, 2])
  end
end
