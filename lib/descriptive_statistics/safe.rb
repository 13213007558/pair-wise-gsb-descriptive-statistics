require "descriptive_statistics/support/convert"
require 'descriptive_statistics/number.rb'
require 'descriptive_statistics/sum.rb'
require 'descriptive_statistics/mean.rb'
require 'descriptive_statistics/median.rb'
require 'descriptive_statistics/mode.rb'
require 'descriptive_statistics/variance.rb'
require 'descriptive_statistics/standard_deviation.rb'
require 'descriptive_statistics/percentile.rb'
require 'descriptive_statistics/percentile_rank.rb'
require 'descriptive_statistics/range.rb'
require 'descriptive_statistics/modified_band_depth.rb'
require 'descriptive_statistics/descriptive_statistics.rb'
require 'descriptive_statistics/stats.rb'
require 'descriptive_statistics/class_methods.rb'

class << DescriptiveStatistics
  def modified_band_depth(curves)
    DescriptiveStatistics::CurveBands.modified_band_depth(curves)
  end

  def modified_band_depth_order(curves)
    DescriptiveStatistics::CurveBands.modified_band_depth_order(curves)
  end

  def sort_by_modified_band_depth(curves)
    DescriptiveStatistics::CurveBands.sort_by_modified_band_depth(
      DescriptiveStatistics::CurveBands.normalize(curves)
    )
  end

  def modified_band_depth_details(curves)
    DescriptiveStatistics::CurveBands.modified_band_depth_details(curves)
  end

  def central_envelope(curves, proportion = 0.5)
    DescriptiveStatistics::CurveBands.central_envelope(curves, proportion)
  end

  alias central_curve_envelope central_envelope
end
