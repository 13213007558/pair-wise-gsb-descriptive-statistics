module DescriptiveStatistics
  def hartigan_dip(
    collection = self,
    simulations: nil,
    random: nil,
    &block
  )
    values = Support::convert(collection, &block)
    return if values.empty?

    result = HartiganDip.result(values)

    unless simulations.nil?
      simulations = simulations.to_int
      raise ArgumentError, "simulations must be a positive integer" if simulations <= 0

      generator = random || Random.new
      exceedances = 0
      simulations.times do
        uniform_sample = Array.new(values.size) { generator.rand }
        simulated_result = HartiganDip.result(uniform_sample)
        simulated_statistic = simulated_result[:statistic]
        exceedances += 1 if simulated_statistic >= result[:statistic]
      end

      p_value = exceedances.fdiv(simulations)
      result[:p_value] = p_value
      result[:simulations] = simulations
      result[:standard_error] = Math.sqrt(p_value * (1.0 - p_value) / simulations)
      result[:random] = generator
    end

    result
  end

  # Pure-Ruby port of Hartigan and Hartigan's AS 217 algorithm.  Indexes in
  # this class are one-based to follow the published algorithm directly.
  class HartiganDip
    def self.result(values)
      statistic, low, high, gcm_indices, lcm_indices = compute(values)
      x = values.map(&:to_f).sort

      {
        statistic: statistic,
        modal_interval: [x[low - 1], x[high - 1]],
        convex_lower_bound_contacts: contact_points(x, gcm_indices, :lower),
        concave_upper_bound_contacts: contact_points(x, lcm_indices, :upper)
      }
    end

    def self.statistic(values)
      dip, = compute(values)
      dip
    end

    def self.compute(values)
      x = values.map(&:to_f).sort
      n = x.size
      dip = 1.0
      low = 1
      high = n
      predecessor = nil
      successor = nil

      if n >= 2 && x[n - 1] != x[0]
        predecessor = build_predecessors(x)
        successor = build_successors(x)
        dip, low, high = find_dip(
          x, predecessor, successor, low, high, dip
        )
      end

      gcm_indices, lcm_indices =
        final_chains(predecessor, successor, low, high, n)
      [dip / (2.0 * n), low, high, gcm_indices, lcm_indices]
    end

    def self.build_predecessors(x)
      n = x.size
      predecessor = Array.new(n + 1)
      predecessor[1] = 1

      (2..n).each do |j|
        predecessor[j] = j - 1
        loop do
          current = predecessor[j]
          previous = predecessor[current]
          break if current == 1 ||
                   (x[j - 1] - x[current - 1]) * (current - previous) <
                     (x[current - 1] - x[previous - 1]) * (j - current)

          predecessor[j] = previous
        end
      end

      predecessor
    end

    def self.build_successors(x)
      n = x.size
      successor = Array.new(n + 1)
      successor[n] = n

      (n - 1).downto(1) do |k|
        successor[k] = k + 1
        loop do
          current = successor[k]
          following = successor[current]
          break if current == n ||
                   (x[k - 1] - x[current - 1]) * (current - following) <
                     (x[current - 1] - x[following - 1]) * (k - current)

          successor[k] = following
        end
      end

      successor
    end

    def self.find_dip(x, predecessor, successor, low, high, initial_dip)
      n = x.size
      dip = initial_dip

      loop do
        gcm = Array.new(n + 1)
        gcm[1] = high
        gcm_position = 1
        while gcm[gcm_position] > low
          gcm_position += 1
          gcm[gcm_position] = predecessor[gcm[gcm_position - 1]]
        end
        gcm_length = gcm_position
        gcm_counter = gcm_length - 1

        lcm = Array.new(n + 1)
        lcm[1] = low
        lcm_position = 1
        while lcm[lcm_position] < high
          lcm_position += 1
          lcm[lcm_position] = successor[lcm[lcm_position - 1]]
        end
        lcm_length = lcm_position
        lcm_counter = 2

        maximum_distance = 0.0
        lower_contact = gcm_counter + 1
        upper_contact = lcm_length

        unless gcm_length == 2 && lcm_length == 2
          loop do
            if gcm[gcm_counter] > lcm[lcm_counter]
              current_gcm = gcm[gcm_counter]
              next_gcm = gcm[gcm_counter + 1]
              current_lcm = lcm[lcm_counter]
              distance = lower_envelope_distance(
                x, current_lcm, next_gcm, current_gcm
              )
              lcm_counter += 1
              if distance >= maximum_distance
                maximum_distance = distance
                lower_contact = gcm_counter + 1
                upper_contact = lcm_counter - 1
              end
            else
              current_gcm = gcm[gcm_counter]
              current_lcm = lcm[lcm_counter]
              previous_lcm = lcm[lcm_counter - 1]
              distance = upper_envelope_distance(
                x, current_gcm, previous_lcm, current_lcm
              )
              gcm_counter -= 1
              if distance >= maximum_distance
                maximum_distance = distance
                lower_contact = gcm_counter + 1
                upper_contact = lcm_counter
              end
            end

            gcm_counter = 1 if gcm_counter < 1
            lcm_counter = lcm_length if lcm_counter > lcm_length
            break if gcm[gcm_counter] == lcm[lcm_counter]
          end
        else
          maximum_distance = 1.0
        end
        break if maximum_distance < initial_dip

        lower_dip, = current_side_dip(
          x, gcm, lower_contact, gcm_length, :lower
        )
        upper_dip, = current_side_dip(
          x, lcm, upper_contact, lcm_length, :upper
        )

        if upper_dip > lower_dip
          new_dip = upper_dip
        else
          new_dip = lower_dip
        end
        dip = new_dip if dip < new_dip

        new_low = gcm[lower_contact]
        new_high = lcm[upper_contact]
        break if low == new_low && high == new_high

        low = new_low
        high = new_high
      end

      [dip, low, high]
    end

    def self.lower_envelope_distance(x, lcm_index, gcm_start, gcm_end)
      numerator = (x[lcm_index - 1] - x[gcm_start - 1]) *
                  (gcm_end - gcm_start)
      denominator = x[gcm_end - 1] - x[gcm_start - 1]
      (lcm_index - gcm_start + 1) - numerator.fdiv(denominator)
    end

    def self.upper_envelope_distance(x, gcm_index, lcm_start, lcm_end)
      numerator = (x[gcm_index - 1] - x[lcm_start - 1]) *
                  (lcm_end - lcm_start)
      denominator = x[lcm_end - 1] - x[lcm_start - 1]
      numerator.fdiv(denominator) -
        (gcm_index - lcm_start - 1).to_f
    end

    def self.current_side_dip(x, chain, start, length, side)
      maximum = 1.0
      index = -1

      (start...length).each do |position|
        chain_start, chain_end = [chain[position], chain[position + 1]].minmax
        next unless chain_end - chain_start > 1

        slope = (chain_end - chain_start).fdiv(
          x[chain_end - 1] - x[chain_start - 1]
        )
        (chain_start..chain_end).each do |data_index|
          value = if side == :lower
            (data_index - chain_start + 1).to_f -
              (x[data_index - 1] - x[chain_start - 1]) * slope
          else
            (x[data_index - 1] - x[chain_start - 1]) * slope -
              (data_index - chain_start - 1).to_f
          end
          next unless maximum < value

          maximum = value
          index = data_index
        end
      end

      [maximum, index]
    end

    def self.final_chains(predecessor, successor, low, high, n)
      return [[], []] if n < 2 || predecessor.nil?

      gcm = []
      index = high
      while index >= low
        gcm << index
        break if index == low
        index = predecessor[index]
      end

      lcm = []
      index = low
      while index <= high
        lcm << index
        break if index == high
        index = successor[index]
      end

      [gcm, lcm]
    end

    def self.contact_points(x, indices, side)
      indices.sort.map do |index|
        probability = side == :lower ? (index - 1).fdiv(x.size) : index.fdiv(x.size)
        [x[index - 1], probability]
      end
    end
  end
end
