require 'safe/spec_helper'

describe "DescriptiveStatistics.hartigan_dip" do
  require 'descriptive_statistics/safe'

  subject { DescriptiveStatistics }

  it "returns nil for an empty collection" do
    expect(subject.hartigan_dip([])).to be_nil
  end

  it "uses the AS 217 minimum for one observation" do
    result = subject.hartigan_dip([7])
    expect(result[:statistic]).to eq(0.5)
    expect(result[:modal_interval]).to eq([7.0, 7.0])
    expect(result[:convex_lower_bound_contacts]).to eq([])
    expect(result[:concave_upper_bound_contacts]).to eq([])
  end

  it "uses 1/(2n) for fewer than four observations" do
    expect(subject.hartigan_dip([7, 9])[:statistic]).to eq(0.25)
    expect(subject.hartigan_dip([7, 7, 9])[:statistic]).to eq(1.0 / 6)
    expect(subject.hartigan_dip([7, 7, 7])[:statistic]).to eq(1.0 / 6)
  end

  it "uses 1/(2n), not zero, for four identical observations" do
    result = subject.hartigan_dip([4, 4, 4, 4])
    expect(result[:statistic]).to eq(0.125)
    expect(result[:modal_interval]).to eq([4.0, 4.0])
    expect(result[:convex_lower_bound_contacts]).to eq([])
    expect(result[:concave_upper_bound_contacts]).to eq([])
  end

  it "computes a fixed four-point uniform reference value" do
    result = subject.hartigan_dip([0, 1, 2, 3])
    expect(result[:statistic]).to eq(0.125)
    expect(result[:modal_interval]).to eq([0.0, 3.0])
    expect(result[:convex_lower_bound_contacts]).to eq([[0.0, 0.0], [3.0, 0.75]])
    expect(result[:concave_upper_bound_contacts]).to eq([[0.0, 0.25], [3.0, 1.0]])
  end

  it "computes a fixed tied bimodal reference value" do
    result = subject.hartigan_dip([1, 1, 2, 2])
    expect(result[:statistic]).to eq(0.25)
    expect(result[:modal_interval]).to eq([2.0, 2.0])
    expect(result[:convex_lower_bound_contacts]).to eq([[2.0, 0.5], [2.0, 0.75]])
    expect(result[:concave_upper_bound_contacts]).to eq([[2.0, 0.75], [2.0, 1.0]])
  end

  it "is invariant to translation and positive proportional rescaling" do
    data = [0, 1, 1, 3, 6, 7]
    transformed = data.map { |value| 2.5 * value + 4.0 }
    original = subject.hartigan_dip(data)
    scaled = subject.hartigan_dip(transformed)

    expect(scaled[:statistic]).to eq(original[:statistic])
    expect(scaled[:modal_interval]).to eq(
      original[:modal_interval].map { |value| 2.5 * value + 4.0 }
    )
  end

  it "mirrors the statistic, interval, and contact coordinates" do
    data = [0, 1, 1, 3, 6, 7]
    center = 8.0
    original = subject.hartigan_dip(data)
    mirrored = subject.hartigan_dip(data.map { |value| center - value })

    mirror_point = proc do |point|
      [center - point.first, 1.0 - point.last]
    end

    expect(mirrored[:statistic]).to be_within(1e-15).of(original[:statistic])
    expect(mirrored[:modal_interval]).to eq(
      original[:modal_interval].reverse.map { |value| center - value }
    )
    mirrored[:convex_lower_bound_contacts].zip(
      original[:concave_upper_bound_contacts].reverse.map(&mirror_point)
    ) do |actual, expected|
      expect(actual.first).to eq(expected.first)
      expect(actual.last).to be_within(1e-15).of(expected.last)
    end
    mirrored[:concave_upper_bound_contacts].zip(
      original[:convex_lower_bound_contacts].reverse.map(&mirror_point)
    ) do |actual, expected|
      expect(actual.first).to eq(expected.first)
      expect(actual.last).to be_within(1e-15).of(expected.last)
    end
  end

  it "is available on enumerable collections" do
    collection = [1, 1, 2, 2].extend(DescriptiveStatistics)
    result = collection.hartigan_dip
    expect(result[:statistic]).to eq(0.25)
  end

  it "computes a deterministic uniform-null Monte Carlo p-value with simulation error" do
    first = subject.hartigan_dip(
      [1, 1, 2, 2], simulations: 100, random: Random.new(12345)
    )
    second = subject.hartigan_dip(
      [1, 1, 2, 2], simulations: 100, random: Random.new(12345)
    )

    expect(first[:p_value]).to eq(0.0)
    expect(first[:simulations]).to eq(100)
    expect(first[:standard_error]).to eq(0.0)
    expect(first[:p_value]).to eq(second[:p_value])
    expect(first[:standard_error]).to eq(second[:standard_error])
  end

  it "rejects zero Monte Carlo simulations" do
    expect { subject.hartigan_dip([1, 2, 3, 4], simulations: 0) }
      .to raise_error(ArgumentError)
  end
end
