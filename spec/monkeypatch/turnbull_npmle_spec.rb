require 'monkeypatch/spec_helper'
require 'descriptive_statistics'

describe DescriptiveStatistics::TurnbullNpmle do
  def direct_npmle(records, dimensions)
    probabilities = Array.new(dimensions, 1.0 / dimensions)
    step = 0.1

    20_000.times do
      gradient = Array.new(dimensions, 0.0)
      records.each do |row|
          value = row.each_with_index.inject(0.0) do |total, (included, index)|
            total + (included ? probabilities[index] : 0.0)
          end
        row.each_with_index do |included, index|
          gradient[index] += 1.0 / value if included
        end
      end

      current_value = likelihood(records, probabilities)
      candidate = nil
      50.times do
        candidate = project_simplex(probabilities.each_with_index.map do |probability, index|
          probability + step * gradient[index]
        end)
        if likelihood(records, candidate) >= current_value
          probabilities = candidate
          step = 0.1
          break
        end
        step *= 0.5
      end
    end

    probabilities
  end

  def project_simplex(values)
    sorted = values.sort.reverse
    cumulative = 0.0
    rho = 0
    theta = -1.0

    sorted.each_with_index do |value, index|
      cumulative += value
      if value - (cumulative - 1.0) / (index + 1).to_f > 0.0
        rho = index
      end
    end
    theta = (sorted[0..rho].inject(0.0, :+) - 1.0) / (rho + 1).to_f

    values.map do |value|
      projected = value - theta
      projected > 0.0 ? projected : 0.0
    end
  end

  def likelihood(records, probabilities)
    records.sum do |row|
      value = row.each_with_index.inject(0.0) do |total, (included, index)|
        total + (included ? probabilities[index] : 0.0)
      end
      Math.log(value)
    end
  end

  it 'estimates a nondegenerate interval-censored Turnbull support' do
    records = [
      {type: :interval, lower: 0, upper: 2, left_closed: false, right_closed: false},
      {type: :interval, lower: 1, upper: 3, left_closed: false, right_closed: false},
      {type: :interval, lower: 2, upper: 4, left_closed: false, right_closed: false}
    ]
    result = described_class.new.fit(records)

    expect(result.converged).to be(true)
    expect(result.support).to eq([
      {lower: 1.0, upper: 2.0, left_closed: false, right_closed: false,
       identifiable: false},
      {lower: 2.0, upper: 3.0, left_closed: false, right_closed: false,
       identifiable: false}
    ])
    expect(result.probabilities).to all(be_between(0.0, 1.0))
    expect(result.probabilities.sum).to be_within(1e-12).of(1.0)
    expect(result.log_likelihood).to be_within(1e-8).of(2.0 * Math.log(0.5))
    expect(result.iterations).to be > 0
    expect(result.active_set_iterations).to be > 0
    expect(result.kkt_residual).to be <= 1e-8
    expect(result.support_size).to eq(5)
    expect(result.active_set_size).to eq(2)

    rows = [
      [true, true, false, false, false],
      [false, true, true, true, false],
      [false, false, false, true, true]
    ]
    reference = direct_npmle(rows, 5)
    expect(result.log_likelihood).to be_within(1e-7).of(likelihood(rows, reference))
    expect(result.probabilities).to all(be_within(2e-6).of(0.5))
  end

  it 'keeps endpoint openness explicit and builds half-open maximal intersections' do
    records = [
      {type: :interval, lower: 0, upper: 1, left_closed: false, right_closed: true},
      {type: :interval, lower: 2, upper: 3, left_closed: false, right_closed: true}
    ]
    result = described_class.new.fit(records)

    expect(result.support).to eq([
      {lower: 0.0, upper: 1.0, left_closed: false, right_closed: true,
       identifiable: false},
      {lower: 2.0, upper: 3.0, left_closed: false, right_closed: true,
       identifiable: false}
    ])
    expect(result.probabilities[0]).to be_within(1e-8).of(0.5)
    expect(result.probabilities[1]).to be_within(1e-8).of(0.5)
    expect(result.support_size).to eq(2)
    expect(result.active_set_size).to eq(2)
  end

  it 'keeps an identifiable shared endpoint separate from non-identifiable intervals' do
    records = [
      {type: :interval, lower: 0, upper: 1, left_closed: false, right_closed: true},
      {type: :interval, lower: 1, upper: 2, left_closed: true, right_closed: false}
    ]
    result = described_class.new.fit(records)

    expect(result.support).to eq([
      {lower: 1.0, upper: 1.0, left_closed: true, right_closed: true,
       identifiable: true}
    ])
    expect(result.probabilities).to eq([1.0])
    expect(result.cdf(1.0)).to eq(1.0)
  end

  it 'reduces exact observations to the empirical distribution' do
    values = [1.0, 2.0, 2.0, 3.0]
    records = values.map { |value| {type: :exact, value: value} }
    result = described_class.new.fit(records)

    expect(result.converged).to be(true)
    expect(result.support).to eq([
      {lower: 1.0, upper: 1.0, left_closed: true, right_closed: true,
       identifiable: true},
      {lower: 2.0, upper: 2.0, left_closed: true, right_closed: true,
       identifiable: true},
      {lower: 3.0, upper: 3.0, left_closed: true, right_closed: true,
       identifiable: true}
    ])
    expect(result.probabilities).to eq([0.25, 0.5, 0.25])
    expect(result.log_likelihood).to be_within(1e-12).of(
      Math.log(0.25) + 2.0 * Math.log(0.5) + Math.log(0.25)
    )
    expect(result.cdf_range(0.0)).to eq(lower: 0.0, upper: 0.0)
    expect(result.cdf(0.0)).to eq(0.0)
    expect(result.cdf(1.0)).to eq(0.25)
    expect(result.cdf(2.0)).to eq(0.75)
    expect(result.cdf(3.0)).to eq(1.0)
  end

  it 'returns sharp CDF bounds around interval mass' do
    records = [
      {type: :left, value: 2.0, right_closed: false},
      {type: :right, value: 2.0, left_closed: false}
    ]
    result = described_class.new.fit(records)

    expect(result.support).to eq([
      {lower: -Float::INFINITY, upper: 2.0, left_closed: false,
       right_closed: false, identifiable: false},
      {lower: 2.0, upper: Float::INFINITY, left_closed: false,
       right_closed: false, identifiable: false}
    ])
    expect(result.probabilities).to all(be_within(1e-10).of(0.5))
    expect(result.cdf_range(2.0)).to eq(lower: 0.5, upper: 0.5)
    expect(result.cdf(2.0)).to eq(0.5)
    expect(result.cdf_range(1.0)).to eq(lower: 0.0, upper: 0.5)
    expect(result.cdf_range(4.0)).to eq(lower: 0.5, upper: 1.0)
  end

  it 'requires endpoint openness explicitly for censored intervals' do
    expect do
      described_class.new.fit([{type: :left, value: 1.0}])
    end.to raise_error(ArgumentError, /right_closed must be declared explicitly/)

    expect do
      described_class.new.fit([
        {type: :interval, lower: -Float::INFINITY, upper: 1.0,
         left_closed: true, right_closed: true}
      ])
    end.to raise_error(ArgumentError, /left_closed cannot be true/)
  end
end
