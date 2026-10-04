require 'safe/spec_helper'

describe "DescriptiveStatistics::ModifiedBandDepth" do
  require 'descriptive_statistics/safe'

  # Crossing curves (not parallel shifts of one another):
  # f0 rises, f1 falls, f2 is constant and always inside their envelope.
  let(:crossing_curves) do
    [
      [0.0, 1.0, 2.0, 3.0],
      [3.0, 2.0, 1.0, 0.0],
      [1.0, 1.0, 1.0, 1.0]
    ]
  end

  subject { DescriptiveStatistics::ModifiedBandDepth.new(crossing_curves) }

  it "ranks the curve staying inside the crossing envelope deepest" do
    expect(subject.depth_values).to eql([0.25, 0.25, 1.0])
    expect(subject.ranking).to eql([2, 0, 1])
  end

  it "exposes the exact numerator and denominator per curve" do
    accounting = subject.depths.map { |d| [d.index, d.inside, d.total] }
    expect(accounting).to eql([[0, 1, 4], [1, 1, 4], [2, 4, 4]])
  end

  it "builds a central envelope from the deepest proportion of curves" do
    envelope = subject.central_envelope(0.5)
    expect(envelope[:indices]).to eql([2, 0])
    expect(envelope[:lower]).to eql([0.0, 1.0, 1.0, 1.0])
    expect(envelope[:upper]).to eql([1.0, 1.0, 2.0, 3.0])
  end

  context "with missing observations" do
    let(:sparse_curves) do
      [
        [1.0, nil, 3.0],
        [0.0, 0.0, 0.0],
        [2.0, 2.0, 2.0],
        [nil, 1.0, 1.0]
      ]
    end

    subject { DescriptiveStatistics::ModifiedBandDepth.new(sparse_curves) }

    it "counts only fully observed (position, curve pair) combinations in the denominator" do
      # For curve 0 the valid combinations are:
      #   pair (1,2) at positions 0 and 2 -> 2 (curve 0 missing at 1)
      #   pair (1,3) at position 2        -> 1 (curve 3 missing at 0, curve 0 at 1)
      #   pair (2,3) at position 2        -> 1
      # Nothing is imputed as zero, so the denominator is 4, not 9.
      depth = subject.depths[0]
      expect(depth.total).to eql(4)
      expect(depth.inside).to eql(1)
      expect(depth.value).to eql(0.25)
    end

    it "never treats missing values as zero" do
      # If nil were read as 0.0, curve 0 would also fall inside the band of
      # pair (2,3) at position 0 and inside pair (1,2) at position 1,
      # scoring 3 instead of 1.
      expect(subject.depths[0].inside).to eql(1)
    end

    it "returns nil envelope positions where no central curve is observed" do
      curves = [
        [0.0, 1.0, 2.0, 3.0],
        [3.0, 2.0, 1.0, 0.0],
        [1.0, nil, 1.0, 1.0]
      ]
      mbd = DescriptiveStatistics::ModifiedBandDepth.new(curves)
      envelope = mbd.central_envelope(0.33)
      expect(envelope[:indices]).to eql([2])
      expect(envelope[:lower]).to eql([1.0, nil, 1.0, 1.0])
      expect(envelope[:upper]).to eql([1.0, nil, 1.0, 1.0])
    end

    it "sorts curves without any valid combination last with a nil depth" do
      curves = crossing_curves + [[nil, nil, nil, nil]]
      mbd = DescriptiveStatistics::ModifiedBandDepth.new(curves)
      expect(mbd.depth_values[3]).to be_nil
      expect(mbd.depths[3].total).to eql(0)
      expect(mbd.ranking.last).to eql(3)
    end
  end

  context "with tied depths" do
    let(:tied_curves) do
      [
        [0.0, 1.0, 2.0, 3.0],
        [3.0, 2.0, 1.0, 0.0],
        [1.0, 1.0, 1.0, 1.0],
        [1.0, 1.0, 1.0, 1.0]
      ]
    end

    subject { DescriptiveStatistics::ModifiedBandDepth.new(tied_curves) }

    it "gives identical curves identical depths" do
      expect(subject.depth_values[2]).to eql(1.0)
      expect(subject.depth_values[3]).to eql(1.0)
    end

    it "keeps the original curve order for ties (stable ranking)" do
      expect(subject.ranking).to eql([2, 3, 0, 1])
    end
  end

  context "validation" do
    it "rejects curves observed on different grids" do
      expect do
        DescriptiveStatistics::ModifiedBandDepth.new([[1.0, 2.0], [1.0, 2.0, 3.0], [1.0, 2.0]])
      end.to raise_error(ArgumentError)
    end

    it "rejects fewer than three curves" do
      expect do
        DescriptiveStatistics::ModifiedBandDepth.new([[1.0], [2.0]])
      end.to raise_error(ArgumentError)
    end

    it "rejects proportions outside (0, 1]" do
      expect { subject.central_envelope(0) }.to raise_error(ArgumentError)
      expect { subject.central_envelope(1.5) }.to raise_error(ArgumentError)
    end
  end
end
