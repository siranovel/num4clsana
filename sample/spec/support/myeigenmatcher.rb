# 固有値・固有ベクトル
module MyEigenMatcher
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
            act_edvec = actual[i][:edvec]
            exp_edval = @expected[i][:edval]
            exp_edvec = @expected[i][:edvec]
            # 固有値
            if (act_edval.round(@n) != exp_edval) then
                ret = false
            end
            # 固有ベクトル
            act_edvec.size.times{|j|
                if act_edvec[j].round(@n) != exp_edvec[j] then
                    ret = false
                end
            }
        end
        return ret
    end
    def failure_message
      "#{@expected} expected but got #{@actual}"
    end
  end
  def is_eigens(expected, n)
    Matcher.new(expected, n)
  end
end

