require 'safe/spec_helper'
require 'descriptive_statistics'

describe DescriptiveStatistics::Turnbull do
  Turnbull = DescriptiveStatistics::Turnbull

  # Independent reference: enumerates every distribution on the elementary
  # sets whose masses are multiples of +step+ and returns the best
  # log likelihood found by direct maximization.
  def brute_force_best_log_likelihood(observations, step)
    endpoints = observations
      .flat_map { |obs| [obs.left, obs.right] }
      .reject { |e| e.infinite? }
      .uniq
      .sort
    sets = [[-Float::INFINITY, endpoints.first]]
    endpoints.each_with_index do |t, j|
      sets << [t, t]
      sets << [t, endpoints[j + 1]] if j + 1 < endpoints.size
    end
    sets << [endpoints.last, Float::INFINITY]

    covered = sets.map do |left, right|
      observations.map do |obs|
        if left == right
          (obs.left < left || (obs.left == left && obs.left_closed?)) &&
            (obs.right > right || (obs.right == right && obs.right_closed?))
        else
          obs.left <= left && obs.right >= right
        end
      end
    end

    best = -Float::INFINITY
    units = (1.0 / step).round
    masses = []
    enumerate = nil
    enumerate = lambda do |slot, remaining|
      if slot == sets.size - 1
        masses << remaining
        log_likelihood = 0.0
        observations.each_with_index do |_obs, obs_index|
          probability = 0.0
          sets.each_index do |i|
            probability += masses[i] * step if covered[i][obs_index]
          end
          log_likelihood = -Float::INFINITY and break if probability <= 0
          log_likelihood += Math.log(probability)
        end
        best = log_likelihood if log_likelihood > best
        masses.pop
      else
        (0..remaining).each do |u|
          masses << u
          enumerate.call(slot + 1, remaining - u)
          masses.pop
        end
      end
    end
    enumerate.call(0, units)
    best
  end

  def masses_by_interval(fit)
    fit.support.map { |entry| [entry.interval, entry.mass] }.to_h
  end

  context "with interval censored observations" do
    # (0,2], (1,3], (2,3]: the likelihood is (a+x)(x+b)(b) with a+x+b = 1,
    # maximized at a = 0, x = b = 1/2, giving log likelihood -2*log(2).
    subject(:fit) do
      Turnbull.new([
        Turnbull.interval(0, 2, left_closed: false, right_closed: true),
        Turnbull.interval(1, 3, left_closed: false, right_closed: true),
        Turnbull.interval(2, 3, left_closed: false, right_closed: true)
      ]).fit
    end

    it "maximizes the likelihood at the analytically known value" do
      expect(fit.log_likelihood).to be_within(1e-6).of(-2 * Math.log(2))
    end

    it "matches a brute-force direct optimization of the likelihood" do
      reference = brute_force_best_log_likelihood(fit.observations, 0.1)
      expect(fit.log_likelihood).to be >= reference - 1e-6
    end

    it "puts half the mass on each maximal overlap interval" do
      masses = masses_by_interval(fit)
      expect(masses.size).to eql(2)
      expect(masses[Turnbull.interval(1, 2, left_closed: false, right_closed: true)]).to be_within(1e-4).of(0.5)
      expect(masses[Turnbull.interval(2, 3, left_closed: false, right_closed: true)]).to be_within(1e-4).of(0.5)
    end

    it "does not collapse unidentifiable mass onto the interval midpoint" do
      fit.support.each do |entry|
        expect(entry.interval.left).to be < entry.interval.right
      end
    end

    it "reports a CDF range inside support intervals and exact values elsewhere" do
      expect(fit.cdf(0.5)).to eql(0.0..0.0)
      expect(fit.cdf(1.5).first).to be_within(1e-9).of(0.0)
      expect(fit.cdf(1.5).last).to be_within(1e-6).of(0.5)
      expect(fit.cdf(2).first).to be_within(1e-6).of(0.5)
      expect(fit.cdf(2).last).to be_within(1e-6).of(0.5)
      expect(fit.cdf(2.5).first).to be_within(1e-6).of(0.5)
      expect(fit.cdf(2.5).last).to be_within(1e-9).of(1.0)
      expect(fit.cdf(3).first).to be_within(1e-9).of(1.0)
      expect(fit.cdf(3).last).to be_within(1e-9).of(1.0)
    end
  end

  context "with a mix of exact and interval censored observations" do
    # exact 1, exact 3, (0,4): likelihood s1*s3*1 with s1+s3 <= 1 is
    # maximized at s1 = s3 = 1/2, giving log likelihood -2*log(2).
    subject(:fit) do
      Turnbull.new([
        Turnbull.exact(1),
        Turnbull.exact(3),
        Turnbull.interval(0, 4, left_closed: false, right_closed: false)
      ]).fit
    end

    it "matches a brute-force direct optimization of the likelihood" do
      reference = brute_force_best_log_likelihood(fit.observations, 0.05)
      expect(fit.log_likelihood).to be >= reference - 1e-6
      expect(fit.log_likelihood).to be_within(1e-6).of(-2 * Math.log(2))
    end

    it "splits the mass evenly between the two exact points" do
      masses = masses_by_interval(fit)
      expect(masses[Turnbull.exact(1)]).to be_within(1e-4).of(0.5)
      expect(masses[Turnbull.exact(3)]).to be_within(1e-4).of(0.5)
    end
  end

  context "with only exact observations" do
    subject(:fit) { Turnbull.new([1, 2, 2, 3].map { |v| Turnbull.exact(v) }).fit }

    it "degenerates to the empirical distribution" do
      masses = masses_by_interval(fit)
      expect(masses.size).to eql(3)
      expect(masses[Turnbull.exact(1)]).to be_within(1e-9).of(0.25)
      expect(masses[Turnbull.exact(2)]).to be_within(1e-9).of(0.5)
      expect(masses[Turnbull.exact(3)]).to be_within(1e-9).of(0.25)
    end

    it "has the empirical log likelihood" do
      expected = Math.log(0.25) + 2 * Math.log(0.5) + Math.log(0.25)
      expect(fit.log_likelihood).to be_within(1e-9).of(expected)
    end

    it "has a point-valued CDF equal to the empirical CDF" do
      expect(fit.cdf(0.5)).to eql(0.0..0.0)
      expect(fit.cdf(1)).to eql(0.25..0.25)
      expect(fit.cdf(1.5)).to eql(0.25..0.25)
      expect(fit.cdf(2)).to eql(0.75..0.75)
      expect(fit.cdf(3)).to eql(1.0..1.0)
      expect(fit.cdf(4)).to eql(1.0..1.0)
    end
  end

  context "with right censored observations" do
    # exact 1, [2,inf), [2,inf): likelihood s*(1-s)^2, maximized at s = 1/3.
    subject(:fit) do
      Turnbull.new([
        Turnbull.exact(1),
        Turnbull.right_censored(2),
        Turnbull.right_censored(2)
      ]).fit
    end

    it "places the remaining mass on the censored ray, not at a point" do
      masses = masses_by_interval(fit)
      expect(masses[Turnbull.exact(1)]).to be_within(1e-4).of(1.0 / 3)
      expect(masses[Turnbull.interval(2, Float::INFINITY, left_closed: true, right_closed: false)])
        .to be_within(1e-4).of(2.0 / 3)
    end

    it "has the analytically known log likelihood" do
      expected = Math.log(1.0 / 3) + 2 * Math.log(2.0 / 3)
      expect(fit.log_likelihood).to be_within(1e-6).of(expected)
    end

    it "leaves the CDF unidentifiable beyond the censoring point" do
      expect(fit.cdf(1).first).to be_within(1e-6).of(1.0 / 3)
      expect(fit.cdf(1).last).to be_within(1e-6).of(1.0 / 3)
      expect(fit.cdf(100).first).to be_within(1e-6).of(1.0 / 3)
      expect(fit.cdf(100).last).to be_within(1e-9).of(1.0)
    end
  end

  context "with a left censored observation" do
    # exact 0, (-inf,-1]: likelihood s*(1-s), maximized at s = 1/2.
    subject(:fit) do
      Turnbull.new([
        Turnbull.exact(0),
        Turnbull.left_censored(-1)
      ]).fit
    end

    it "splits the mass across the censoring bound" do
      masses = masses_by_interval(fit)
      expect(masses[Turnbull.exact(0)]).to be_within(1e-4).of(0.5)
      expect(masses[Turnbull.interval(-Float::INFINITY, -1, left_closed: false, right_closed: true)])
        .to be_within(1e-4).of(0.5)
      expect(fit.log_likelihood).to be_within(1e-6).of(2 * Math.log(0.5))
    end

    it "has an identifiable CDF outside the censored ray" do
      expect(fit.cdf(0).first).to be_within(1e-9).of(1.0)
      expect(fit.cdf(0).last).to be_within(1e-9).of(1.0)
      expect(fit.cdf(-0.5).first).to be_within(1e-6).of(0.5)
      expect(fit.cdf(-0.5).last).to be_within(1e-6).of(0.5)
      expect(fit.cdf(-1.5).first).to be_within(1e-9).of(0.0)
      expect(fit.cdf(-1.5).last).to be_within(1e-6).of(0.5)
    end
  end

  context "with explicitly open endpoints" do
    # (0,2) open on both sides plus an exact observation at 2: the open
    # interval does not cover 2, so the likelihood is s*(1-s), maximized
    # at s = 1/2.
    subject(:fit) do
      Turnbull.new([
        Turnbull.interval(0, 2, left_closed: false, right_closed: false),
        Turnbull.exact(2)
      ]).fit
    end

    it "keeps the open interval separate from the exact endpoint" do
      masses = masses_by_interval(fit)
      expect(masses.size).to eql(2)
      expect(masses[Turnbull.interval(0, 2, left_closed: false, right_closed: false)])
        .to be_within(1e-4).of(0.5)
      expect(masses[Turnbull.exact(2)]).to be_within(1e-4).of(0.5)
    end

    it "respects openness in the CDF" do
      expect(fit.cdf(0)).to eql(0.0..0.0)
      expect(fit.cdf(1).first).to be_within(1e-9).of(0.0)
      expect(fit.cdf(1).last).to be_within(1e-6).of(0.5)
      expect(fit.cdf(2).first).to be_within(1e-9).of(1.0)
      expect(fit.cdf(2).last).to be_within(1e-9).of(1.0)
    end
  end

  context "diagnostics" do
    subject(:fit) do
      Turnbull.new([
        Turnbull.interval(0, 2, left_closed: false, right_closed: true),
        Turnbull.interval(1, 3, left_closed: false, right_closed: true),
        Turnbull.interval(2, 3, left_closed: false, right_closed: true)
      ], tolerance: 1e-10).fit
    end

    it "reports convergence details" do
      diagnostics = fit.diagnostics
      expect(diagnostics[:converged]).to be true
      expect(diagnostics[:iterations]).to be > 0
      expect(diagnostics[:iterations]).to be <= diagnostics[:max_iterations]
      expect(diagnostics[:max_probability_change]).to be < diagnostics[:tolerance]
      expect(diagnostics[:log_likelihood]).to eql(fit.log_likelihood)
    end

    it "has a non-decreasing log likelihood trace" do
      trace = fit.diagnostics[:log_likelihood_trace]
      expect(trace.size).to eql(fit.diagnostics[:iterations])
      trace.each_cons(2) do |previous, current|
        expect(current).to be >= previous - 1e-12
      end
    end
  end

  context "input handling" do
    it "accepts plain numerics as exact observations" do
      fit = Turnbull.new([1, 2, 2, 3]).fit
      expect(masses_by_interval(fit)[Turnbull.exact(2)]).to be_within(1e-9).of(0.5)
    end

    it "fits lazily when queried" do
      expect(Turnbull.new([1, 2, 2, 3]).cdf(2)).to eql(0.75..0.75)
    end

    it "rejects empty data" do
      expect { Turnbull.new([]) }.to raise_error(ArgumentError)
    end

    it "rejects reversed endpoints" do
      expect { Turnbull.interval(3, 1) }.to raise_error(ArgumentError)
    end

    it "rejects open degenerate intervals" do
      expect { Turnbull.interval(2, 2, right_closed: false) }.to raise_error(ArgumentError)
    end

    it "rejects closed infinite endpoints" do
      expect { Turnbull.interval(1, Float::INFINITY, right_closed: true) }.to raise_error(ArgumentError)
      expect { Turnbull.interval(-Float::INFINITY, 1, left_closed: true) }.to raise_error(ArgumentError)
    end
  end
end
