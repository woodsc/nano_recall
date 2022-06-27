#represents a sample, generally from a fastq file.
require "#{File.dirname(__FILE__)}/settings.rb"

class Sample
  attr_accessor :id, :source_filepath, :genes

  def initialize()
    @id = nil
    @source_filepath = nil
    @genes = []
  end

  #produces stats on the sample
  def print_stats()
    puts "Statistics #{@id}"
    @genes.each do |gene|
      gene_seq_cnt = gene.clades.map(){|c| c.alignments.size()}.sum()
      puts "  Gene #{gene.name}:  #{gene.clades.size()} clades, #{gene_seq_cnt} sequences."
      gene.clades.sort(){|a,b| b.alignments.size() <=> a.alignments.size() }.each do |clade|
        puts "    Clade #{clade.id}:  #{clade.alignments.size()} sequences."
      end
    end
  end


end


#sequences from a gene, organized by clade.
class Gene
  attr_accessor :name, :clades, :aa_coverage, :hxb2_standard, :trim_range

  def initialize(name: , hxb2_standard: )
    @name = name
    @clades = []
    @hxb2_standard = hxb2_standard
    @aa_coverage = []
    0.upto((@hxb2_standard.nucleotides.size() / 3) - 1) do |i|
      @aa_coverage[i] = 0
    end
    @trim_range = (0 .. @hxb2_standard.nucleotides.size() - 1)
  end


  #returns the clades organized by subtype instead of exact reference.
  def subtypes()
    data = []
    @clades.map(){|clade| clade.subtype}.uniq().each do |subtype|
      clades = @clades.find_all() {|clade| clade.subtype == subtype}
      new_clade = Clade.new()
      new_clade.id = subtype
      new_clade.subtype = subtype
      new_clade.ref_standard = clades[0].ref_standard
      new_clade.hxb2_standard = clades[0].hxb2_standard
      clades.each do |clade|
        new_clade.alignments += clade.alignments
      end
      data << new_clade
    end

    return data
  end

end

class Clade
  attr_accessor :id, :subtype, :alignments, :ref_standard, :hxb2_standard

  def initialize()
    @id = nil #Ref.A1.AU.03.PS1044_Day0.DQ676872
    @subtype = nil #A, B, C, AE, etc
    @alignments = []
    @ref_standard = nil
    @hxb2_standard = nil
  end


  #Find the most common variants.  We could potentially use this for data correction.
  #Also slightly biases the first few sequences, but very slightly.
  def find_variants(filename) #experimental
    variant_list = [] #[seq, cnt]

    @alignments.each() do |a|
      no_match = true
      variant_list.each_with_index do |variant, v_i|
        a_seq = a.seq_clean()
        diff = 0
        #ignore inserts and deletions.
        0.upto(variant[0].size() - 1) do |i|
          if(variant[0][i] != a_seq[i] and variant[0][i] != '-' and a_seq[i] != '-')
            diff += 1
            #puts "Diff #{variant[0][i]} #{a_seq[i]}"
          end
        end
        puts "Variant #{v_i} and #{a.seq.id} differ by #{diff}"
        if(diff == 0)
          no_match = false
          variant[1] += 1
          #Also, we can fill in any missing pieces of the variant.
          0.upto(variant.size() - 1) do |i|
            if(variant[0][i] == '-' and a_seq[i] != '-')
              variant[0][i] = a_seq[i]
              puts "fill"
            end
          end
          break
        end
      end
      if(no_match)
        variant_list << [a.seq_clean(), 1]
      end

    end

    variant_list.sort(){|a,b| b[1] <=> a[1]}.each do |variant|
      puts variant[1]
    end

    File.open(filename, 'w') do |file|
      file.puts "Sequence,Count"
      variant_list.sort(){|a,b| b[1] <=> a[1]}.each do |variant|
        file.puts variant.join(',')
      end
    end

  end



end
