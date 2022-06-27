require "#{File.dirname(__FILE__)}/../../lib/alignment/_alignment.rb"

data = []

def match_score(align)
  score = 0.0
  0.upto(align[0].size() - 1) do |i|
    if(align[0][i] == '$')
      score += 0.0
    elsif(align[0][i] == align[1][i])
      score += 1.0
    elsif(align[0][i] == '-' or align[1][i] == '-')
      score += 0.5
    end
  end

  score -= [align[0].count('-'), align[1].count('-')].min() * 3
  #score -= 6  if(align[0].include?('-') and align[1].include?('-'))
  return score
end

File.open("testdata.csv", 'r') do |file|
  file.each_line do |line|
    data << line.strip.split(',')
  end
end

criteria = [
  [0,0,0,[]],
  [0,1,0,[]],
  [0,2,0,[]],
  [0,3,0,[]],
  [0,4,0,[]],
  [1,0,0,[]],
  [1,1,0,[]],
  [1,2,0,[]],
  [1,3,0,[]],
  [1,4,0,[]],
  [2,0,0,[]],
  [2,1,0,[]],
  [2,2,0,[]],
  [2,3,0,[]],
  [2,4,0,[]],
  [3,0,0,[]],
  [3,1,0,[]],
  [3,2,0,[]],
  [3,3,0,[]],
  [3,4,0,[]],
  [4,0,0,[]],
  [4,1,0,[]],
  [4,2,0,[]],
  [4,3,0,[]],
  [4,4,0,[]],
]


#Okay, we need to test different alignment parameters and find the best one.
data.each do |dat|
  criteria.each do |crit|
    if(dat[1] and dat[2])
      alignment = align_it('$' + dat[2] + '$', '$' + dat[1] + '$', crit[0], crit[1])
      crit[3] << alignment
      crit[2] += match_score(alignment)
      #puts match_score(alignment)
    end
  end
end

criteria.each do |crit|
  puts "#{crit[0]} #{crit[1]}:  #{crit[2]}"
end
puts "Best:"
criteria.sort!(){|a,b| b[2] <=> a[2]}
criteria[0 .. 4].each() do |crit|
  puts "#{crit[0]} #{crit[1]}:  #{crit[2]}"
end

STDIN.gets

criteria[0][3].each do |alignment|
  puts alignment
  puts match_score(alignment)
end
