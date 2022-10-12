require "#{File.dirname(__FILE__)}/../settings.rb"
require "#{File.dirname(__FILE__)}/../sample.rb"


#experimental, disabled by default
def reject_poor_aligned_ends(alignment: )
  good = true

  if(Settings['reject-poor-aligned-ends'])
    window = Settings['reject-poor-aligned-ends-window']
    size_threshold = Settings['reject-poor-aligned-ends-size-threshold']
    match_threshold = Settings['reject-poor-aligned-ends-match-threshold']

    rejected = 0

    alignment.trim_start.upto(alignment.trim_end) do |i|
      mat, tot = 0, 0
      (i - window).upto(alignment.trim_end) do |j|
        next if(j < alignment.trim_start or j > alignment.trim_end)
        next if(alignment.std.nucleotides[j] == '-')
        tot += 1
        mat += 1 if(match_nuc(alignment.std.nucleotides[j], alignment.seq.nucleotides[j]))
        break if(tot == window)
      end
      match_perc = mat.to_f / tot.to_f

      if(match_perc < match_threshold)
        rejected += 1
      else
        break #done
      end
    end

    alignment.trim_end.downto(alignment.trim_start) do |i|
      mat, tot = 0, 0
      (i + window).downto(alignment.trim_start) do |j|
        next if(j < alignment.trim_start or j > alignment.trim_end)
        next if(alignment.std.nucleotides[j] == '-')
        tot += 1
        mat += 1 if(match_nuc(alignment.std.nucleotides[j], alignment.seq.nucleotides[j]))
        break if(tot == window)
      end
      match_perc = mat.to_f / tot.to_f

      if(match_perc < match_threshold)
        rejected += 1
      else
        break #done
      end
      #puts "#{i} #{alignment.std.nucleotides[i]}  #{match_perc}  #{match_perc < 0.6 ? 'BAD' : ''}" if(alignment.seq.id == '063514a7-041c-4ba7-b7a2-9551cea48dbf')
      #puts "#{i} #{alignment.std.nucleotides[i]}  #{match_perc}  #{match_perc < 0.6 ? 'BAD' : ''}"
    end


    #clean this up later.
    if(rejected > 0)
      #reject
      seq_size = alignment.trim_end - alignment.trim_start
      bad_perc = rejected.to_f / seq_size.to_f
      #puts "#{alignment.seq.id}:  #{(bad_perc * 100).round()}% poor edge alignment"
      if(bad_perc >= size_threshold)
        #puts "#{alignment.seq.id}:  REJECTED"
        good = false #rejected
      else
        #puts "I GUESS ITS FINE"
      end
    end

  end

  return good
end
