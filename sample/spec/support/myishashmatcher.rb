# ハッシュ
module MyIsHashMatcher
  class Matcher
    def initialize(expected, n)
      @expected = expected
      @n = n
    end
    def matches?(actual)
        @actual = actual
        ret = (@actual.size == @expected.size) ? true : false
        @expected.each do |k, v|
            if @expected[k] != @actual[k].round(@n) then
                ret = false
            end
        end
        return ret
    end
    def failure_message
      "#{@expected} expected but got #{@actual}"
    end
  end
  def is_hash(expected, n)
    Matcher.new(expected, n)
  end
end

