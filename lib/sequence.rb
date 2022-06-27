=begin
sequence.rb's goals are to encapsulate a single sequence in a memory efficient
manner.

Well, first lets just get it providing useful methods.  We can worry about
memory usage when/if it becomes a problem.

We'll need some way to keep alignments represented?  IE, quality scores should
be aligned as well.

Do we want a seperate alignment class maybe?  Or chain things?  Hrm..
Lets keep it simple for now.  We can add complexity when we need it.

=end

require "#{File.dirname(__FILE__)}/settings.rb"

class Sequence
  attr_accessor :quality, :nucleotides, :annotations, :id, :source_filepath
  attr_accessor :reverse_complemented, :reference

  def initialize(nuc: nil, qual: nil, id: nil, reference: false, source_filepath: nil)
    @quality = (qual ? qual : [])
    @nucleotides = (nuc ? nuc : "")
    @id = (id ? id : nil)
    @annotations = {}
    @reference = reference
    @source_filepath = source_filepath
    @reverse_complemented = false
  end

  #returns a new reverse complement sequence.
  def reverse_complement()
    rc = self.copy()
    rc.reverse_complemented = !rc.reverse_complemented
    rc.quality = rc.quality.reverse()
    rc.nucleotides = rc.nucleotides.reverse()
    rc.nucleotides = rc.nucleotides.tr('ATGC','0123').tr('0123','TACG')
    return rc
  end

  def copy()
    seq = Sequence.new()
    seq.id = @id
    seq.reference = @reference
    seq.source_filepath = @source_filepath
    seq.annotations = @annotations
    seq.reverse_complemented = @reverse_complemented
    seq.quality = @quality
    seq.nucleotides = @nucleotides
    return seq
  end

  def to_s()
    str = "#<Sequence:#{self.object_id}>\n"
    str += @nucleotides + "\n"
    str += "  id: #{@id}\n"
    str += "  reference seq?: #{@reference}  reverse_complemented?:  #{@reverse_complemented}\n"
    str += "  source_filepath: #{@source_filepath}\n"
    str += "  annotations: #{@annotations}\n"
    str += "  length: #{@nucleotides.size()}  has_qual?:  #{(@quality and @quality.size() > 0)}\n"

    return str
  end

  def compare_match_cnt(seq) #unused?
    matches = 0
    minsize = [@nucleotides.size(), seq.nucleotides.size()].min()
    0.upto(minsize - 1) do |i|
      #if(@nucleotides[i] == seq.nucleotides[i])
      if(match_nuc(@nucleotides[i], seq.nucleotides[i]))
        matches += 1
      end
    end
    return matches
  end

  #Only counts gaps in the original sequence.
  def compare_match_perc(seq) #unused?
    matches = 0
    minsize = [@nucleotides.size(), seq.nucleotides.size()].min()
    gaps = 0
    0.upto(minsize - 1) do |i|
      if(@nucleotides[i] == '-')
        gaps += 1
      #elsif(@nucleotides[i] == seq.nucleotides[i])
      elsif(match_nuc(@nucleotides[i], seq.nucleotides[i]))
        matches += 1
      end
    end
    return matches.to_f / (minsize - gaps).to_f
  end





end
