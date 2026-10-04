Descriptive Statistics
===============
[![Gem Version](https://badge.fury.io/rb/descriptive_statistics.svg)](http://badge.fury.io/rb/descriptive_statistics)
[![Build Status](https://travis-ci.org/thirtysixthspan/descriptive_statistics.svg?branch=master)](https://travis-ci.org/thirtysixthspan/descriptive_statistics)

Overview
--------

This gem adds methods to the Enumerable module to allow easy calculation of basic
descriptive statistics of Numeric sample data in collections that have included Enumerable such as Array, Hash, Set, and Range. The statistics that can be calculated are:
* Number
* Sum
* Mean
* Median
* Mode
* Variance
* Standard Deviation
* Percentile
* Percentile Rank
* Descriptive Statistics
* Quartiles
* Distance Covariance / Correlation


When requiring DescriptiveStatistics, the Enumerable module is monkey patched so
that the statistical methods are available on any instance of a class that has included Enumerable. For example with an Array:
```ruby
require 'descriptive_statistics'
data = [2,6,9,3,5,1,8,3,6,9,2]
# => [2, 6, 9, 3, 5, 1, 8, 3, 6, 9, 2]
data.number
# => 11.0
data.sum
# => 54.0
data.mean
# => 4.909090909090909
data.median
# => 5.0
data.variance
# => 7.7190082644628095
data.standard_deviation
# => 2.778310325442932
data.percentile(30)
# => 3.0
data.percentile(70)
# => 6.0
data.percentile_rank(8)
# => 81.81818181818183
data.mode
# => 2
data.range
# => 8
data.descriptive_statistics
# => {:number=>11.0,
  :sum=>54,
  :variance=>7.7190082644628095,
  :standard_deviation=>2.778310325442932,
  :min=>1,
  :max=>9,
  :mean=>4.909090909090909,
  :mode=>2,
  :median=>5.0,
  :range=>8.0,
  :q1=>2.5,
  :q2=>5.0,
  :q3=>7.0}
```

and with other types of objects:
```ruby
require 'set'
require 'descriptive_statistics'
{:a=>1, :b=>2, :c=>3, :d=>4, :e=>5}.mean #Hash
# => 3.0
Set.new([1,2,3,4,5]).mean #Set
# => 3.0
(1..5).mean #Range
# => 3.0
```

including instances of your own classes, when an `each` method is provided that
creates an Enumerator over the sample data to be operated upon:
```ruby
class Foo
  include Enumerable
  attr_accessor :bar, :baz, :bat

  def each
    [@bar, @baz, @bat].each
  end
end

foo = Foo.new
foo.bar = 1
foo.baz = 2
foo.bat = 3
foo.mean
# => 2.0

```

or:
```ruby
class Foo
  include Enumerable
  attr_accessor :bar, :baz, :bat

  def each
    Enumerator.new do |y|
      y << @bar
      y << @baz
      y << @bat
    end
  end
end

foo = Foo.new
foo.bar = 1
foo.baz = 2
foo.bat = 3
foo.mean
# => 2.0
```

and even Structs:
```ruby
Scores = Struct.new(:sally, :john, :peter)
bowling = Scores.new
bowling.sally = 203
bowling.john = 134
bowling.peter = 233
bowling.mean
# => 190.0
```

All methods optionally take blocks that operate on object values. For example:
```ruby
require 'descriptive_statistics'
LineItem = Struct.new(:price, :quantity)
cart = [ LineItem.new(2.50, 2), LineItem.new(5.10, 9), LineItem.new(4.00, 5) ]
total_items = cart.sum(&:quantity)
# => 16
total_price = cart.sum{ |i| i.price * i.quantity }
# => 70.9
```

Note that you can extend DescriptiveStatistics on individual objects by
requiring DescriptiveStatistics safely, thus avoiding the monkey patch. For example:
```ruby
require 'descriptive_statistics/safe'
data = [2,6,9,3,5,1,8,3,6,9,2]
# => [2, 6, 9, 3, 5, 1, 8, 3, 6, 9, 2]
data.extend(DescriptiveStatistics)
# => [2, 6, 9, 3, 5, 1, 8, 3, 6, 9, 2]
data.number
# => 11.0
data.sum
# => 54
```

Or, if you prefer leaving your collection pristine, you can create a
Stats object that references your collection:
```ruby
require 'descriptive_statistics/safe'
data = [1, 2, 3, 4, 5, 1]
# => [1, 2, 3, 4, 5, 1]
stats = DescriptiveStatistics::Stats.new(data)
# => [1, 2, 3, 4, 5, 1]
stats.class
# => DescriptiveStatistics::Stats
stats.mean
# => 2.6666666666666665
stats.median
# => 2.5
stats.mode
# => 1
data << 2
# => [1, 2, 3, 4, 5, 1, 2]
data << 2
# => [1, 2, 3, 4, 5, 1, 2, 2]
stats.mode
# => 2
```

Or you call the statistical methods directly:
```ruby
require 'descriptive_statistics/safe'
# => true
DescriptiveStatistics.mean([1,2,3,4,5])
# => 3.0
DescriptiveStatistics.mode([1,2,3,4,5])
# => 1
DescriptiveStatistics.variance([1,2,3,4,5])
# => 2.0
```

Or you can use [Refinements](http://www.ruby-doc.org/core/doc/syntax/refinements_rdoc.html) (available in Ruby >= 2.1) to augment any class that mixes in the Enumerable module. Refinements are lexically scoped and so the statistical methods will only be available in the file where they are used. Note that the lexical scope can be limited to a Class or Module, but only applies to code in that file. This approach provides a great deal of protection against introducing conflicting modifications to Enumerable while retaining the convenience of the monkey patch approach.
```ruby
require 'descriptive_statistics/refinement'

class SomeServiceClass
  using DescriptiveStatistics::Refinement.new(Array)

  def self.calculate_something(array)
    array.standard_deviation
  end
end

[1,2,3].standard_deviation
# => NoMethodError: undefined method `standard_deviation' for [1, 2, 3]:Array

SomeServiceClass.calculate_something([1,2,3])
#=> 0.816496580927726
```


Ruby on Rails
-------------

To use DescriptiveStatistics with Ruby on Rails add DescriptiveStatistics to your Gemfile, requiring the safe extension.
```ruby
source 'https://rubygems.org'

gem 'rails', '4.1.7'
gem 'descriptive_statistics', '~> 2.4.0', :require => 'descriptive_statistics/safe'
```
Then after a `bundle install`, you can extend DescriptiveStatistics on an individual collection and call the statistical methods as needed.
```ruby
users = User.all.extend(DescriptiveStatistics)
mean_age = users.mean(&:age) # => 19.428571428571427
mean_age_in_dog_years = users.mean { |user| user.age / 7.0 } # => 2.7755102040816326
```

Distance Dependence
-------------------
Distance covariance and distance correlation measure dependence between two
random vectors, even when the two sides have different dimensions or Pearson
correlation is close to zero. Each input row represents one paired observation:
the first element is the X vector and the second element is the Y vector. Plain
Numeric values are treated as one-dimensional vectors.
```ruby
require 'descriptive_statistics/safe'

pairs = [
  [[0, 0], [0, 0, 0]],
  [[1, 0], [1, 0, 0]],
  [[2, 0], [2, 0, 0]],
  [[3, 0], [3, 0, 0]]
]

DescriptiveStatistics.distance_covariance(pairs)
# => 0.9013878188659973
DescriptiveStatistics.distance_correlation(pairs)
# => 1.0

# The same call is available through the monkey patch or refinements:
pairs.extend(DescriptiveStatistics).distance_correlation

# Two zipped collections and blocks are also supported:
DescriptiveStatistics.distance_covariance([0, 1, 2, 3], [0, 1, 2, 3])
records = [{ :x => 0, :y => 0 }, { :x => 1, :y => 1 }]
DescriptiveStatistics.distance_covariance(records) { |record| [record[:x], record[:y]] }
```

Use `distance_dependence_statistics` to get the ordinary sample estimates and
the U-centered squared estimates together. The ordinary values are the usual
nonnegative distance covariance and correlation, normalized by n squared. The
`u_centered_squared_*` values use the U-centered inner-product normalization
`n * (n - 3)`; the covariance and variance components are unbiased squared
estimates. The squared correlation ratio is bias corrected rather than a
guaranteed nonnegative quantity, so negative finite-sample values are returned
as calculated and are never clipped to zero.
```ruby
DescriptiveStatistics.distance_dependence_statistics(pairs)
# => {
#      :distance_covariance => ...,
#      :distance_correlation => ...,
#      :squared_distance_covariance => ...,
#      :squared_distance_correlation => ...,
#      :u_centered_squared_distance_covariance => ...,
#      :u_centered_squared_distance_correlation => ...,
#      ...
#    }

DescriptiveStatistics.distance_covariance(pairs, :estimator => :u_centered_squared)
DescriptiveStatistics.distance_correlation(pairs, :estimator => :u_centered_squared)
```

By default, `missing_rows: :complete_rows` drops a paired row when either side
contains `nil` or `Float::NAN`, and reports the number dropped as
`:dropped_rows`. Use `missing_rows: :raise` to reject missing rows instead.
Non-finite values other than NaN, non-Numeric values, empty vectors, and
within-side dimension changes always raise `ArgumentError`. The ordinary
estimator requires at least two complete rows; the U-centered estimator
requires at least four.

For a paired permutation test, complete Y vectors are permuted while X rows
remain fixed. The default statistic is `:u_centered_squared`; pass
`estimator: :sample` to use the ordinary (square-root scale) distance
covariance. The Monte Carlo p value includes the observed ordering as
`(exceedances + 1) / (repetitions + 1)`.
Provide `seed:` for reproducibility, or any object responding to `rand(maximum)`
as `random:`.
```ruby
DescriptiveStatistics.distance_covariance_test(
  pairs,
  :repetitions => 999,
  :seed => 1234
)
# => { :statistic => ..., :p_value => ..., :repetitions => 999, ... }
```

Notes
-----
* All methods return a Float object except for `mode`, which will return a Numeric object from the collection. `mode` will always return nil for empty collections.
* All methods return nil when the collection is empty, except for `number`, which returns 0.0. This is a different behavior than [ActiveSupport's Enumerable monkey patch of sum](http://apidock.com/rails/Enumerable/sum), which by deafult returns the Fixnum 0 for empty collections. You can change this behavior by specifying the default value returned for empty collections all at once:
  ```ruby
require 'descriptive_statistics'
[].mean
# => nil
[].sum
# => nil
DescriptiveStatistics.empty_collection_default_value = 0.0
# => 0.0
[].mean
# => 0.0
[].sum
# => 0.0
  ```
  or one at a time:
  ```ruby
require 'descriptive_statistics'
[].mean
# => nil
[].sum
# => nil
DescriptiveStatistics.sum_empty_collection_default_value = 0.0
# => 0.0
[].mean
# => nil
[].sum
# => 0.0
DescriptiveStatistics.mean_empty_collection_default_value = 0.0
# => 0.0
[].mean
# => 0.0
[].sum
# => 0.0
  ```

* The scope of this gem covers [Descriptive Statistics](http://en.wikipedia.org/wiki/Descriptive_statistics) and not Inferential Statistics. From wikipedia:

  > Descriptive statistics is the discipline of quantitatively describing the main features of a collection of information, or the quantitative description itself. Descriptive statistics are distinguished from inferential statistics, in that descriptive statistics aim to summarize a sample, rather than use the data to learn about the population that the sample of data is thought to represent.

  Thus, all statistics calculated herein describe only the values in the collection. Where this makes a practical difference is in the calculation of variance (and thus the standard deviation which is derived from variance). We use the [population variance](http://en.wikipedia.org/wiki/Variance#Population_variance) to calculate the variance of the values in the collection. If the values in your collection represent a sampling from a larger population of values, and your goal is to estimate the population variance from your sample, you should use the inferential statistic, [sample variance](http://en.wikipedia.org/wiki/Variance#Sample_variance). However, the calculation of the sample variance is outside the scope of this gem's functionality.


Ports
-----
[Javascript](http://github.com/FGRibreau/descriptive_statistics)

[Python](http://github.com/gleicon/py_descriptive_statistics)

[Go](https://github.com/gleicon/go-descriptive-statistics)

[Elixir](https://github.com/pusewicz/descriptive_statistics)

[Clojure](https://github.com/nickmcdonnough/descriptivestatistics)

License
-------
Copyright (c) 2010-2014
Derrick Parkhurst (derrick.parkhurst@gmail.com),
Gregory Brown (gregory.t.brown@gmail.com),
Daniel Farrell (danielfarrell76@gmail.com),
Graham Malmgren,
Guy Shechter,
Charlie Egan (charlieegan3@googlemail.com)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
