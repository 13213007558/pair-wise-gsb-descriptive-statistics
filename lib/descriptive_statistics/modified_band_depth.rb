module DescriptiveStatistics
  module CurveBands
    module_function

    def modified_band_depth(curves)
      modified_band_depth_details(curves).map { |details| details.fetch(:depth) }
    end

    def modified_band_depth_order(curves)
      depths = modified_band_depth(curves)

      depths.each_with_index.sort_by { |depth, index|
        depth.nil? ? [1, 0.0, index] : [0, -depth, index]
      }.map(&:last)
    end

    def sort_by_modified_band_depth(curves)
      modified_band_depth_order(curves).map { |index| curves[index] }
    end

    def modified_band_depth_details(curves)
      curves = normalize(curves)
      validate_common_grid!(curves)

      curves.map.with_index do |curve, curve_index|
        other_indices = (0...curves.size).to_a - [curve_index]
        bands = other_indices.combination(2).map do |first_index, second_index|
          observed_positions = []
          contained_positions = []

          curve.each_with_index do |value, point|
            first_value = curves[first_index][point]
            second_value = curves[second_index][point]
            next unless observed_value?(value) &&
                        observed_value?(first_value) &&
                        observed_value?(second_value)

            observed_positions << point
            lower = [first_value, second_value].min
            upper = [first_value, second_value].max
            if value >= lower && value <= upper
              contained_positions << point
            end
          end

          {
            :curve_indices => [first_index, second_index],
            :observed_positions => observed_positions,
            :contained_positions => contained_positions,
            :depth => band_depth(observed_positions, contained_positions)
          }
        end

        included_bands = bands.reject { |band| band.fetch(:observed_positions).empty? }
        depth = if included_bands.empty?
                  nil
                else
                  included_bands.sum { |band| band.fetch(:depth) } / included_bands.size
                end

        {
          :curve_index => curve_index,
          :depth => depth,
          :observed_positions => included_bands.sum { |band|
            band.fetch(:observed_positions).size
          },
          :contained_positions => included_bands.sum { |band|
            band.fetch(:contained_positions).size
          },
          :bands => bands
        }
      end
    end

    def central_envelope(curves, proportion)
      curves = normalize(curves)
      validate_common_grid!(curves)
      validate_proportion!(proportion)
      selected_curves = central_envelope_curves(curves, proportion)
      point_count = curves.first&.size || 0

      lower = Array.new(point_count)
      upper = Array.new(point_count)

      point_count.times do |point|
        values = selected_curves.map { |curve| curve[point] }.select {
          |value| observed_value?(value)
        }
        next if values.empty?

        lower[point] = values.min
        upper[point] = values.max
      end

      [lower, upper]
    end

    def normalize(curves)
      Array(curves).map { |curve| Array(curve) }
    end

    def central_envelope_curves(curves, proportion)
      depths = modified_band_depth(curves)
      defined_depths = depths.reject(&:nil?)
      return [] if curves.empty? || proportion.zero?

      if defined_depths.empty?
        raise ArgumentError, "no curves have an observed band"
      end

      selected_count = (defined_depths.size * proportion).ceil
      selected_count = 1 if selected_count.zero?
      selected_indices = modified_band_depth_order(curves).first(selected_count)

      selected_indices.map { |index| curves[index] }
    end

    def validate_proportion!(proportion)
      return if proportion.is_a?(Numeric) && proportion >= 0 && proportion <= 1

      raise ArgumentError, "proportion must be between 0 and 1"
    end

    def validate_common_grid!(curves)
      point_counts = curves.map(&:size).uniq
      return if point_counts.size <= 1

      raise ArgumentError, "all curves must have the same number of grid positions"
    end

    def observed_value?(value)
      !value.nil? && !(value.is_a?(Float) && value.nan?)
    end

    def band_depth(observed_positions, contained_positions)
      return nil if observed_positions.empty?

      contained_positions.size.to_f / observed_positions.size
    end
  end

  def modified_band_depth(curves = nil)
    CurveBands.modified_band_depth(curves.nil? ? self : curves)
  end

  def modified_band_depth_order(curves = nil)
    CurveBands.modified_band_depth_order(curves.nil? ? self : curves)
  end

  def sort_by_modified_band_depth(curves = nil)
    CurveBands.sort_by_modified_band_depth(
      CurveBands.normalize(curves.nil? ? self : curves)
    )
  end

  def modified_band_depth_details(curves = nil)
    CurveBands.modified_band_depth_details(curves.nil? ? self : curves)
  end

  def central_envelope(proportion: 0.5)
    unless proportion.is_a?(Numeric) && proportion >= 0 && proportion <= 1
      raise ArgumentError, "proportion must be between 0 and 1"
    end

    CurveBands.central_envelope(self, proportion)
  end

  alias_method :central_curve_envelope, :central_envelope
end
