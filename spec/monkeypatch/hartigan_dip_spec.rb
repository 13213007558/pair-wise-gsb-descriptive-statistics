require 'monkeypatch/spec_helper'

describe "Array#hartigan_dip" do
  require 'descriptive_statistics'

  it "returns the Hartigan dip result" do
    result = [1, 1, 2, 2].hartigan_dip
    expect(result[:statistic]).to eq(0.25)
    expect(result[:modal_interval]).to eq([2.0, 2.0])
  end
end
