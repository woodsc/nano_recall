#builds and represents an alignment.
require "#{File.dirname(__FILE__)}/sequence.rb"
require "#{File.dirname(__FILE__)}/alignment/_alignment"
require "#{File.dirname(__FILE__)}/settings.rb"

class Alignment
  attr_accessor :seq, :std, :insertions, :deletions, :trim_start, :trim_end, :rev_comp, :censors
  attr_accessor :cache_seq_clean #set to nil to reset cache

  def initialize()
    @seq = nil
    @std = nil
    @insertions = [] #collection of ranges.  MUST BE SORTED.
    @deletions = [] #collection of ranges.  MUST BE SORTED.
    @trim_start = nil #index of sequence start (first nucleotide).
    @trim_end = nil   #index of sequence end (last nucleotide)
    @cache_seq_clean = nil
    @rev_comp = false
    @censors = []  #collection of ranges to censor  (based on raw alignment)
  end

  #aligns and stores the alignment
  def align(seq:, std:, gap_init: Settings['align-gap-init'], gap_extend: Settings['align-gap-extend'])
    raw = align_it(std.nucleotides, seq.nucleotides, gap_init, gap_extend)

    @rev_comp = true if(seq.reverse_complemented)

    @std = std.copy()
    @seq = seq.copy()

    @std.nucleotides = raw[0]
    @seq.nucleotides = raw[1]

    self.recalc_indels()
  end

  def recalc_indels()
    @insertions = []
    @deletions = []
    @trim_start = nil
    @trim_end = nil
    ins = nil
    del = nil
    #scan for trim locations, deletions, insertions.
    0.upto(@std.nucleotides.size() - 1) do |i|
      std_nuc = @std.nucleotides[i]
      seq_nuc = @seq.nucleotides[i]
      if(@trim_start == nil and std_nuc != '-')
        @trim_start = i
        del = i if(del.nil? and seq_nuc == '-')
      elsif(@trim_start != nil and std_nuc == '-' and seq_nuc != '-')
        #insertion
        ins = i if(ins.nil?)
      elsif(@trim_start != nil and std_nuc != '-' and seq_nuc == '-')
        #deletion
        del = i if(del.nil?)
      end

      if(!ins.nil? and std_nuc != '-')
        @insertions.push((ins .. i - 1))
        ins = nil
      end
      if(!del.nil? and seq_nuc != '-')
        @deletions.push((del .. i - 1))
        del = nil
      end
    end

    #oh, ran to the end, I guess this insertion is actually the trim location.
    if(ins)
      @trim_end = ins - 1
    else
      @trim_end = @seq.nucleotides.size() - 1
    end
    #Eh, a deletion I guess?
    if(del)
      @deletions.push((del .. @seq.nucleotides.size() - 1))
    end
  end

  #returns a sequence trimmed && insertions removed.
  def seq_clean(recalc=false)
    if(recalc or @cache_seq_clean.nil?)
      str = @seq.nucleotides + ''
      #apply censors
      @censors.each do |rng|
        rng.each do |i|
          str[i] = '-'
        end
      end

      str = str[0 .. @trim_end]

      @insertions.reverse.each do |ins|
        next if(ins.first >= str.size())
        str[ins] = ''
      end

      str = str[@trim_start .. -1]
      @cache_seq_clean = str
    end
    return @cache_seq_clean
  end

  #returns a sequence trimmed && insertions removed.
  def std_clean()
    str = @std.nucleotides[0 .. @trim_end]

    @insertions.reverse.each do |ins|
      next if(ins.first >= str.size())
      str[ins] = ''
    end

    str = str[@trim_start .. -1]
    return str
  end

  #outputs sequence nicely, trimmed && insertions removed.
  def details()
    str = ""
    str += ">#{@std.id}\n"
    str += self.std_clean() + "\n"
    str += ">#{@seq.id}\n"
    str += self.seq_clean() + "\n"

    str += "Insertions: " + @insertions.map(){|e|
      if(e.size() == 1)
        self.loc_ref(e.first).to_s + ':' + @seq.nucleotides[e]
      else
        self.loc_ref(e.first).to_s + '-' + (self.loc_ref(e.first) + e.size()).to_s + ':' + @seq.nucleotides[e]
      end
    }.join(' ') + "\n"

    str += "Deletions: " + @deletions.map(){|e|
      if(e.size() == 1)
        self.loc_ref(e.first).to_s
      else
        self.loc_ref(e.first).to_s + '-' + (self.loc_ref(e.first) + e.size()).to_s
      end
    }.join(' ') + "\n"
  end

  #takes a raw index pos, and returns the pos relative to standard.
  def loc_ref(i)
    pos = i
    @insertions.each do |ins|
      if(ins.first < i) #has to be the original value i we are comparing.
        pos -= ins.size()
      else
        break #no need to continue since insertions are ordered
      end
    end
    return pos - @trim_start
  end

  #takes a standard index pos and returns the pos in the raw sequence.
  def std_ref(i)
    pos = i + @trim_start

    @insertions.each do |ins|
      if(ins.first < pos)
        pos += ins.size()
      else
        break #no need to continue since insertions are ordered
      end
    end

    return pos
  end

  #returns the quality value at a position (absolute).
  #We want this to be fast...
  def qual_at(pos: nil, insert: nil)
    if(insert.nil?)
      modpos = pos
      #alter position according to deletions
      @deletions.each do |del|
        if(del.first <= pos and del.last >= pos)
          return nil
        elsif(del.last < pos)
          modpos -= del.size()
        end
      end

      return @seq.quality[modpos]
    end
  end
=begin
  #Good test code.
  std = Sequence.new(nuc: "CCCATGCATGGCATGCCCATGGGGCATGCGGG", )
  seq = Sequence.new(nuc: "CCCTTTATGCATGCATGCATGCATGCTTTGGG", qual: [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31])
  alignment = Alignment.new()
  alignment.align(std: std, seq: seq, gap_init: 1)
  puts alignment.details
  0.upto(alignment.std.nucleotides.size() - 1) do |i|
    puts "#{i}:\t#{alignment.std.nucleotides[i]} #{alignment.seq.nucleotides[i]}  #{alignment.qual_at(pos: i)}"
  end
  exit()
=end



  def match_perc()
    matches = 0
    minsize = [@std.nucleotides.size(), @seq.nucleotides.size()].min()
    gaps = 0
    0.upto(minsize - 1) do |i|
      if(@std.nucleotides[i] == '-')
        gaps += 1
      #elsif(@std.nucleotides[i] == @seq.nucleotides[i])
      elsif(match_nuc(@std.nucleotides[i], @seq.nucleotides[i]))
        matches += 1
      end
    end
    return matches.to_f / (minsize - gaps).to_f
  end

  def match_perc_gaps_okay() #unused?
    matches = 0
    minsize = [@std.nucleotides.size(), @seq.nucleotides.size()].min()
    gaps = 0
    0.upto(minsize - 1) do |i|
      if(@std.nucleotides[i] == '-' or @seq.nucleotides[i] == '-')
        gaps += 1
      elsif(match_nuc(@std.nucleotides[i], @seq.nucleotides[i]))
        matches += 1
      end
    end
    return matches.to_f / (minsize - gaps).to_f
  end

  def to_s()
    puts ">#{@std.id}"
    puts @std.nucleotides
    puts ">#{@seq.id}"
    puts @seq.nucleotides
  end

end
