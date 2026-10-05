require 'descriptive_statistics/safe'

class DeterministicRandom
  def initialize(values)
    @values = values
  end

  def rand
    raise 'random source exhausted' if @values.empty?
    @values.shift
  end
end

describe "DescriptiveStatistics hartigan dip" do
  let(:bimodal) { [1.0, 1.0, 2.0, 2.0] }
  let(:sample) { [0.0, 2.0, 3.0, 5.0, 6.0] }
  let(:reflected_sample) { [0.0, 1.0, 3.0, 4.0, 6.0] }
  let(:statfaculty) do
    [30,33,35,36,37,37,39,39,39,39,39,40,40,40,40,41,42,43,43,43,
     44,44,45,45,46,46,47,47,48,48,48,49,50,50,51,52,52,53,53,53,
     53,53,54,54,57,57,59,60,60,60,61,61,61,61,62,62,62,62,63,66,
     70,72,72]
  end

  it "computes the independent reference statistic and contacts" do
    result = DescriptiveStatistics.hartigan_dip(statfaculty)

    expect(result[:dip]).to be_within(1e-14).of(0.05952380952380952)
    expect(result[:modal_interval]).to eq([39.0, 54.0])
    expect(result[:modal_interval_indices]).to eq([7, 44])
    expect(result[:gcm_indices]).to eq([44, 43, 38, 7])
    expect(result[:lcm_indices]).to eq([7, 11, 15, 42, 44])
    expect(result[:gcm_points]).to eq([54.0, 54.0, 53.0, 39.0])
    expect(result[:lcm_points]).to eq([39.0, 39.0, 40.0, 53.0, 54.0])
    expect(result[:gcm_contact_points]).to eq([[54.0, 43.0 / 63],
                                              [54.0, 42.0 / 63],
                                              [53.0, 37.0 / 63],
                                              [39.0, 6.0 / 63]])
    expect(result[:lcm_contact_points]).to eq([[39.0, 7.0 / 63],
                                              [39.0, 11.0 / 63],
                                              [40.0, 15.0 / 63],
                                              [53.0, 42.0 / 63],
                                              [54.0, 44.0 / 63]])
  end

  it "returns the maximal dip for two repeated points" do
    result = DescriptiveStatistics.hartigan_dip(bimodal)

    expect(result[:dip]).to eq(0.25)
    expect(result[:modal_interval]).to eq([2.0, 2.0])
    expect(result[:modal_interval_indices]).to eq([3, 4])
    expect(result[:gcm_indices]).to eq([4, 3])
    expect(result[:lcm_indices]).to eq([3, 4])
  end

  it "defines minimum values for small and all-equal samples" do
    expect(DescriptiveStatistics.hartigan_dip([7.0])[:dip]).to eq(0.5)
    expect(DescriptiveStatistics.hartigan_dip([7.0, 7.0])[:dip]).to eq(0.25)
    expect(DescriptiveStatistics.hartigan_dip([7.0, 7.0, 7.0])[:dip]).to eq(1.0 / 6.0)
    expect(DescriptiveStatistics.hartigan_dip([7.0] * 8)[:dip]).to eq(1.0 / 16.0)
    expect(DescriptiveStatistics.hartigan_dip([])).to be_nil

    expect(DescriptiveStatistics.hartigan_dip([7.0], true)[:dip]).to eq(0.0)
    expect(DescriptiveStatistics.hartigan_dip([7.0, 7.0], true)[:dip]).to eq(0.0)
    expect(DescriptiveStatistics.hartigan_dip([7.0, 7.0, 7.0], true)[:dip]).to eq(0.0)
    expect(DescriptiveStatistics.hartigan_dip([7.0] * 8, true)[:dip]).to eq(0.0)
  end

  it "preserves dip under translation and positive rescaling" do
    transformed = sample.map { |value| 3.0 * value + 4.0 }
    result = DescriptiveStatistics.hartigan_dip(sample)
    transformed_result = DescriptiveStatistics.hartigan_dip(transformed)

    expect(transformed_result[:dip]).to eq(result[:dip])
    expect(transformed_result[:modal_interval]).to eq([19.0, 22.0])
    expect(transformed_result[:modal_interval_indices]).to eq(result[:modal_interval_indices])
    expect(transformed_result[:gcm_indices]).to eq(result[:gcm_indices])
    expect(transformed_result[:lcm_indices]).to eq(result[:lcm_indices])
  end

  it "preserves dip under reflection and reverses the modal interval" do
    result = DescriptiveStatistics.hartigan_dip(sample)
    reflected_result = DescriptiveStatistics.hartigan_dip(reflected_sample)

    expect(reflected_result[:dip]).to be_within(1e-15).of(result[:dip])
    expect(reflected_result[:modal_interval]).to eq([0.0, 1.0])
    expect(reflected_result[:modal_interval_indices]).to eq([1, 2])
    expect(reflected_result[:gcm_indices]).to eq([2, 1])
    expect(reflected_result[:lcm_indices]).to eq([1, 2])
  end

  it "sorts unsorted input" do
    result = DescriptiveStatistics.hartigan_dip(sample.reverse)

    expect(result[:dip]).to eq(DescriptiveStatistics.hartigan_dip(sample)[:dip])
    expect(result[:modal_interval_indices]).to eq([4, 5])
  end

  it "computes a reproducible Monte Carlo p-value with caller-owned randomness" do
    random_values = [
      0.05,0.10,0.50,0.90,
      0.20,0.30,0.40,0.80,
      0.10,0.20,0.70,0.75
    ]
    result = DescriptiveStatistics.hartigan_dip(
      [0.1,0.2,0.35,0.9], false, 3,
      DeterministicRandom.new(random_values.dup)
    )

    expect(result[:p_value]).to eq(1.0)
    expect(result[:simulations]).to eq(3)
    expect(result[:simulation_standard_error]).to eq(0.0)

    zero_result = DescriptiveStatistics.hartigan_dip(
      bimodal, false, 3,
      DeterministicRandom.new(random_values.dup)
    )
    expect(zero_result[:p_value]).to eq(0.0)
    expect(zero_result[:simulation_standard_error]).to eq(0.0)
  end

  it "reports binomial Monte Carlo standard error" do
    random_values = [
      0.05,0.10,0.50,0.90,
      0.20,0.30,0.40,0.80,
      0.10,0.20,0.70,0.75,
      0.01,0.02,0.98,0.99
    ]
    result = DescriptiveStatistics.hartigan_dip(
      [0.1,0.2,0.8,0.9], false, 4,
      DeterministicRandom.new(random_values)
    )

    expect(result[:p_value]).to eq(0.25)
    expect(result[:simulation_standard_error]).to be_within(1e-14).of(0.21650635094610965)
  end

  it "validates Monte Carlo inputs" do
    expect {
      DescriptiveStatistics.hartigan_dip(sample, false, 0, Random.new)
    }.to raise_error(ArgumentError, /positive integer/)
    expect {
      DescriptiveStatistics.hartigan_dip(sample, false, 10, Object.new)
    }.to raise_error(ArgumentError, /respond to rand/)
  end

  it "does not run simulations unless requested" do
    result = DescriptiveStatistics.hartigan_dip(sample)

    expect(result[:p_value]).to be_nil
    expect(result[:simulations]).to eq(0)
    expect(result[:simulation_standard_error]).to be_nil
  end
end
