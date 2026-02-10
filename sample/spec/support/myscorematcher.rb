# 主成分得点
module MyScoreMatcher
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
            act_score = actual[i][:score]
            exp_edval = @expected[i][:edval]
            exp_score = @expected[i][:score]
            # 固有値
            if (act_edval.round(@n) != exp_edval) then
                ret = false
            end
            # 主成分得点
            act_score.size.times{|j|
                if act_score[j].round(@n) != exp_score[j] then
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
  class Matcher2
    def initialize(expected, n)
      @expected = expected
      @n = n
    end
    def matches?(actual)
        ret = true
        @actual = actual
        act_g1 = actual[:G1]
        act_g2 = actual[:G2]
        # グループ1
        act_g1.size.times{|i|
            if act_g1[i].round(@n) != @expected[:G1][i] then
                ret = false
            end
        }
        # グループ2
        act_g2.size.times{|i|
            if act_g2[i].round(@n) != @expected[:G2][i] then
                ret = false
            end
        }
        return ret
    end
    def failure_message
      "#{@expected} expected but got #{@actual}"
    end
  end
  def is_score(expected, n)
    Matcher.new(expected, n)
  end
  def is_score2(expected, n)
    Matcher2.new(expected, n)
  end
end

