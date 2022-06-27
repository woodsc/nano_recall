require "#{File.dirname(__FILE__)}/../settings.rb"
require "#{File.dirname(__FILE__)}/../sample.rb"

#Not sure how much this really helps.
def trim_bad_ends(alignment: )
  if(Settings['trim-bad-ends'])
    clip_limit = 999999
    if(alignment.trim_end.nil? or alignment.trim_start.nil?)
      clip_limit = (alignment.clean_seq().size() / 4).to_i
    else
      clip_limit = ((alignment.trim_end - alignment.trim_start) / 4).to_i
    end

    clip_bases = 0
    trim = 0
    if(!alignment.trim_end.nil?())
      trim = (alignment.std.nucleotides.size() - alignment.trim_end) - 1
    end

    while(clip_bases < clip_limit)
      clip = 60 + clip_bases
      tmp_align = Alignment.new()
      tmp_align.std = Sequence.new(nuc: alignment.std.nucleotides[-(trim + clip), 60])
      tmp_align.seq = Sequence.new(nuc: alignment.seq.nucleotides[-(trim + clip), 60])

      if(tmp_align.match_perc_gaps_okay() < 0.50)
        #puts "Bad end (#{trim}, #{clip_bases} of #{clip_limit}):  #{alignment.match_perc()} #{tmp_align.match_perc()}"
        alignment.seq.nucleotides[-(trim + (clip - 30)), 30] = ('-' * 30)

        clip_bases += 30
      else
        break
      end
    end

    if(clip_bases >= clip_limit)
      alignment.cache_seq_clean = nil
      alignment.recalc_indels()
      return false #This sequence seems borked, lets not use it.
    elsif(clip_bases > 0)
      if(trim > 0) #get rid of anything past the trim range.
        alignment.seq.nucleotides[-trim, trim] = ''
        alignment.std.nucleotides[-trim, trim] = ''
      end

      #clean sequence end
      (alignment.std.nucleotides.size() - 1).downto(0) do |i|
        if(alignment.std.nucleotides[i] == '-' and alignment.seq.nucleotides[i] == '-')
          alignment.std.nucleotides[i] = ''
          alignment.seq.nucleotides[i] = ''
        end
      end

      alignment.cache_seq_clean = nil
      alignment.recalc_indels()
    end


    #probably should also trim the start, but for now ignore...

  end
  return true
end
