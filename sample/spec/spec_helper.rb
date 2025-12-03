require 'rspec'

Dir["./spec/support/**/*.rb"].each do |f|
    require f
end
RSpec.configure do |config|
  config.include MyEigenMatcher
  config.include MyContriMatcher
  config.include MyScoreMatcher
end

