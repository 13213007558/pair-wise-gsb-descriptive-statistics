module DescriptiveStatistics

  def hartigan_dip(collection = self, min_is_zero = false,
                   simulations = nil, random = nil, &block)
    values = Support.convert(collection, &block).sort
    return DescriptiveStatistics.hartigan_dip_empty_collection_default_value if values.empty?

    result = HartiganDipAlgorithm.result(values, min_is_zero)

    result = {
      :dip => result[:dip],
      :n => values.size,
      :modal_interval => [values[result[:low]], values[result[:high]]],
      :modal_interval_indices => [result[:low] + 1, result[:high] + 1],
      :gcm_indices => result[:gcm].map { |index| index + 1 },
      :lcm_indices => result[:lcm].map { |index| index + 1 },
      :gcm_points => result[:gcm].map { |index| values[index] },
      :lcm_points => result[:lcm].map { |index| values[index] },
      :gcm_contact_points => result[:gcm].map { |index|
        [values[index], index.to_f / values.size]
      },
      :lcm_contact_points => result[:lcm].map { |index|
        [values[index], (index + 1).to_f / values.size]
      },
      :min_is_zero => min_is_zero,
      :null_distribution => "uniform",
      :p_value => nil,
      :simulations => 0,
      :simulation_standard_error => nil
    }
    HartiganDipAlgorithm.add_uniform_p_value(result, values, simulations,
                                             random, min_is_zero)
  end
  alias_method :dip, :hartigan_dip unless method_defined?(:dip)

  def hartigan_dip_test(collection = self, simulations = nil, random = nil,
                        min_is_zero = false, &block)
    values = Support.convert(collection, &block).sort
    return DescriptiveStatistics.hartigan_dip_test_empty_collection_default_value if values.empty?

    hartigan_dip(values, min_is_zero, simulations, random)
  end
  alias_method :dip_test, :hartigan_dip_test unless method_defined?(:dip_test)

  module HartiganDipAlgorithm
    extend self

    def add_uniform_p_value(result, values, simulations, random, min_is_zero)
      return result if simulations.nil?
      unless simulations.is_a?(Integer) && simulations > 0
        raise ArgumentError, "simulations must be a positive integer"
      end
      unless random.respond_to?(:rand)
        raise ArgumentError, "random must respond to rand"
      end

      n = values.size
      if n <= 3
        result[:p_value] = 1.0
        result[:simulations] = simulations
        result[:simulation_standard_error] = 0.0
        return result
      end

      threshold = result[:dip]
      exceedances = 0
      simulations.times do
        sample = Array.new(n) { uniform_order_statistic(random, n) }.sort!
        null_result = result(sample, min_is_zero)
        exceedances += 1 if null_result[:dip] >= threshold
      end

      p_value = exceedances.to_f / simulations
      result[:p_value] = p_value
      result[:simulations] = simulations
      result[:simulation_standard_error] = Math.sqrt(p_value * (1.0 - p_value) / simulations)
      result
    end

    def result(values, min_is_zero)
      hartigan_dip_result(values, min_is_zero)
    end

    private

    def uniform_order_statistic(random, n)
      random.rand.to_f
    end

    def hartigan_dip_result(values, min_is_zero)
      n = values.size
      low = 0
      high = n - 1
      dip = min_is_zero ? 0.0 : 1.0

      return terminal_dip_result(values, low, high, dip.to_f / (2.0 * n)) if n < 2

      return terminal_dip_result(values, low, high, dip.to_f / (2.0 * n)) if values.last == values.first

      previous_gcm = build_previous_gcm(values)
      next_lcm = build_next_lcm(values)
      gcm = []
      lcm = []
      ig = nil
      ih = nil

      loop do
        gcm = [high]
        gcm << previous_gcm[gcm.last] while gcm.last > low
        length_gcm = gcm.size
        ix = length_gcm - 2

        lcm = [low]
        lcm << next_lcm[lcm.last] while lcm.last < high
        length_lcm = lcm.size
        iv = 1

        ig = length_gcm - 1
        ih = length_lcm - 1

        distance = 0.0

        unless length_gcm == 2 && length_lcm == 2
          loop do
            if gcm[ix] > lcm[iv]
              gcm_next = gcm[ix + 1]
              distance_here =
                (lcm[iv] - gcm_next + 1).to_f -
                ((values[lcm[iv]] - values[gcm_next]) *
                 (gcm[ix] - gcm_next).to_f /
                 (values[gcm[ix]] - values[gcm_next]))
              iv += 1
              if distance_here >= distance
                distance = distance_here
                ig = ix + 1
                ih = iv - 1
              end
            else
              lcm_previous = lcm[iv - 1]
              distance_here =
                ((values[gcm[ix]] - values[lcm_previous]) *
                 (lcm[iv] - lcm_previous).to_f /
                 (values[lcm[iv]] - values[lcm_previous])) -
                (gcm[ix] - lcm_previous - 1).to_f
              ix -= 1
              if distance_here >= distance
                distance = distance_here
                ig = ix + 1
                ih = iv
              end
            end

            ix = 0 if ix < 0
            iv = length_lcm - 1 if iv > length_lcm - 1
            break if gcm[ix] == lcm[iv]
          end
        else
          distance = min_is_zero ? 0.0 : 1.0
        end

        break if distance < dip
        lower_dip = 0.0
        (ig...length_gcm - 1).each do |j|
          segment_start = gcm[j + 1]
          segment_end = gcm[j]
          next unless segment_end - segment_start > 1 &&
                      values[segment_end] != values[segment_start]

          slope = (segment_end - segment_start).to_f /
                  (values[segment_end] - values[segment_start])
          (segment_start..segment_end).each do |index|
            candidate =
              (index - segment_start + 1).to_f -
              (values[index] - values[segment_start]) * slope
            lower_dip = [lower_dip, candidate].max
          end
        end

        upper_dip = 0.0
        (ih...length_lcm - 1).each do |j|
          segment_start = lcm[j]
          segment_end = lcm[j + 1]
          next unless segment_end - segment_start > 1 &&
                      values[segment_end] != values[segment_start]

          slope = (segment_end - segment_start).to_f /
                  (values[segment_end] - values[segment_start])
          (segment_start..segment_end).each do |index|
            candidate =
              (values[index] - values[segment_start]) * slope -
              (index - segment_start - 1).to_f
            upper_dip = [upper_dip, candidate].max
          end
        end

        dipnew = [lower_dip, upper_dip].max
        dip = dipnew if dipnew > dip

        break if low == gcm[ig] && high == lcm[ih]

        low = gcm[ig]
        high = lcm[ih]
      end

      {
        :dip => dip / (2.0 * n),
        :low => low,
        :high => high,
        :gcm => gcm,
        :lcm => lcm
      }
    end

    def terminal_dip_result(values, low, high, dip)
      contacts = [low, high].uniq
      {
        :dip => dip,
        :low => low,
        :high => high,
        :gcm => contacts,
        :lcm => contacts
      }
    end

    def build_previous_gcm(values)
      previous = Array.new(values.size)
      previous[0] = 0
      (1...values.size).each do |j|
        previous[j] = j - 1
        loop do
          middle = previous[j]
          before = previous[middle]
          break if middle == 0 ||
                    (values[j] - values[middle]) * (middle - before) <
                    (values[middle] - values[before]) * (j - middle)

          previous[j] = before
        end
      end
      previous
    end

    def build_next_lcm(values)
      following = Array.new(values.size)
      following[values.size - 1] = values.size - 1
      (values.size - 2).downto(0) do |k|
        following[k] = k + 1
        loop do
          middle = following[k]
          after = following[middle]
          break if middle == values.size - 1 ||
                    (values[k] - values[middle]) * (middle - after) <
                    (values[middle] - values[after]) * (k - middle)

          following[k] = after
        end
      end
      following
    end
  end

end
