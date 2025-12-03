# 寄与率
module MyContriMatcher
  class Matcher
    def initialize(expected, n)
      @expected = expected
      @n = n
    end
    def matches?(actual)
        ret = true
        @actual = actual
        actual.size.times do |i|
            act_edval = actual[i][:edval]
            act_cr = actual[i][:cr]
            act_ccr = actual[i][:ccr]
            # 固有値
            if (act_edval.round(@n) != @expected[i][:edval]) then
                ret = false
            end
            # 寄与率
            if (act_cr.round(@n) != @expected[i][:cr]) then
                ret = false
            end
            # 累積寄与率
            if (act_ccr.round(@n) != @expected[i][:ccr]) then
                ret = false
            end
        end
        return ret
    end
    def failure_message
      "#{@expected} expected but got #{@actual}"
    end
  end
  def is_contri(expected, n)
    Matcher.new(expected, n)
  end
end

