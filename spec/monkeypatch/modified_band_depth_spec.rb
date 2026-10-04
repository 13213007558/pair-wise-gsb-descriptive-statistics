require 'monkeypatch/spec_helper'

describe "DescriptiveStatistics modified band depth" do
  require 'descriptive_statistics'

  it "orders crossing curves by containment in other curve envelopes" do
    curves = [
      [-1, 4, -1],
      [4, -1, 4],
      [1, 1, 1]
    ]

    expect(curves.modified_band_depth).to eq([0.0, 0.0, 1.0])
    expect(curves.modified_band_depth_order).to eq([2, 0, 1])
    expect(curves.sort_by_modified_band_depth).to eq([
      curves[2], curves[0], curves[1]
    ])
  end

  it "preserves input order for tied depths" do
    curves = [
      [1, 1, 1],
      [1, 1, 1],
      [1, 1, 1]
    ]

    expect(curves.modified_band_depth).to all(eq(1.0))
    expect(curves.modified_band_depth_order).to eq([0, 1, 2])
  end

  it "only counts grid positions where all three curves are observed" do
    curves = [
      [1, nil, 1, nil],
      [1, 0, 1, 0],
      [1, 2, 1, 2],
      [0, 1, 2, 3],
      [nil, 0, nil, 0]
    ]

    details = curves.modified_band_depth_details

    expect(details[0].fetch(:depth)).to eq(1.0)
    expect(details[0].fetch(:bands).map { |band| band.fetch(:observed_positions) }).to eq([
      [0, 2], [0, 2], [], [0, 2], [], []
    ])
    band_without_positions = details[0].fetch(:bands).find {
      |band| band.fetch(:curve_indices) == [2, 4]
    }
    expect(band_without_positions.fetch(:observed_positions)).to eq([])
    expect(band_without_positions.fetch(:depth)).to be_nil
    first_band = details[1].fetch(:bands).find {
      |band| band.fetch(:curve_indices) == [0, 2]
    }
    expect(first_band.fetch(:observed_positions)).to eq([0, 2])
    expect(first_band.fetch(:depth)).to eq(1.0)
    partial_band = details[1].fetch(:bands).find {
      |band| band.fetch(:curve_indices) == [3, 4]
    }
    expect(partial_band.fetch(:observed_positions)).to eq([1, 3])
    expect(partial_band.fetch(:depth)).to eq(1.0)
  end

  it "returns nil depth for curves without any usable three-curve positions" do
    curves = [
      [nil, nil, 1],
      [1, nil, nil],
      [nil, 1, nil],
      [1, 1, 1]
    ]

    expect(curves.modified_band_depth).to eq([nil, nil, nil, nil])
    expect(curves.modified_band_depth_order).to eq([0, 1, 2, 3])
  end

  it "builds the central envelope from a stable proportion of deepest curves" do
    curves = [
      [-1, 4, -1],
      [4, -1, 4],
      [1, 1, 1],
      [0, 0, 0]
    ]

    expect(curves.central_envelope).to eq([[0, 0, 0], [1, 1, 1]])
  end

  it "does not zero-fill missing envelope positions" do
    curves = [
      [2, nil, 2],
      [2, 2, nil],
      [nil, 2, 2],
      [0, 0, 0]
    ]

    expect(curves.central_envelope).to eq([
      [2, 2, 2],
      [2, 2, 2]
    ])
  end

  it "rejects curves on different grids" do
    expect { [[1], [1, 2]].modified_band_depth }.to raise_error(
      ArgumentError, "all curves must have the same number of grid positions"
    )
  end

  it "supports small curve collections" do
    expect([].modified_band_depth).to eq([])
    expect([[1, 2]].modified_band_depth).to eq([nil])
    expect([[1, 2], [2, 3]].modified_band_depth).to eq([nil, nil])
    expect([].central_envelope).to eq([[], []])
  end

  it "builds an empty envelope for proportion zero" do
    curves = [[-1, 4, -1], [4, -1, 4], [1, 1, 1]]

    expect(curves.central_envelope(proportion: 0)).to eq(
      [[nil, nil, nil], [nil, nil, nil]]
    )
  end

  it "rejects envelope proportions outside 0 through 1" do
    curves = [[-1, 4, -1], [4, -1, 4], [1, 1, 1]]

    expect { curves.central_envelope(proportion: 1.5) }.to raise_error(
      ArgumentError, "proportion must be between 0 and 1"
    )
  end
end
