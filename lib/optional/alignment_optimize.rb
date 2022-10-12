require "#{File.dirname(__FILE__)}/../settings.rb"
require "#{File.dirname(__FILE__)}/../sample.rb"

=begin
1) For each alignment, search for an indel.
2) Realign to most common (and second most common?) sequence in the region.
  Replace if its a better match.
=end

def alignment_optimize(gene: )
  culm_time_a = 0.0
  culm_time_b = 0.0
  culm_time_c = 0.0
  culm_time_d = 0.0
  culm_time_e = 0.0
  culm_time_f = 0.0
  culm_time_g = 0.0
  culm_time_h = 0.0
  culm_time_i = 0.0

  log = []
  padsize = Settings['alignment-optimize-pad-size']
  if(Settings['alignment-optimize'])
    indels = []
    alignments = []
    gene.clades.each do |clade|
      clade.alignments.each do |alignment|
        alignments << alignment
      end
    end


    alignments.each_with_index do |alignment, a_i|
      indels = (alignment.insertions + alignment.deletions).sort(){|a,b| a.first <=> b.first}

      #indels must be processed in reverse order, as each replacement could modify the indexes.
      indels.reverse().each do |indel|

        ts = Time.now
        #Get alignment region.
        frag_rng = (indel.first - 6 .. indel.last + 6)
        while(alignment.std.nucleotides[frag_rng][0] == '-' and frag_rng.first > 0)
          frag_rng = (frag_rng.first - 1 .. frag_rng.last)
        end
        while(alignment.std.nucleotides[frag_rng][-1] == '-' and frag_rng.last < alignment.std.nucleotides.size())
          frag_rng = (frag_rng.first .. frag_rng.last + 1)
        end
        clean_rng = (alignment.loc_ref(frag_rng.first) .. alignment.loc_ref(frag_rng.last))

        next if(clean_rng.first < 0)
        next if(frag_rng.first == 0 or frag_rng.last >= alignment.std.nucleotides.size())

        frag_std = alignment.std.nucleotides[frag_rng]
        frag_seq = alignment.seq.nucleotides[frag_rng]

        culm_time_g += Time.now - ts

        #Find two most common sequences in the alignment region.
        clean_freq = {}
        clean_freq.default = 0

        #slowest
        ts = Time.now
        alignments.each do |sub_alignment|
          st = sub_alignment.std_ref(clean_rng.first)
          ed = sub_alignment.std_ref(clean_rng.last)

          substd = sub_alignment.std.nucleotides[st .. ed]
          subseq = sub_alignment.seq.nucleotides[st .. ed]

          if(substd and subseq and !substd.include?('-') and !subseq.include?('-'))
            clean_freq[subseq] += 1
          end
        end
        culm_time_a += Time.now - ts

        ts = Time.now
        common_subseqs = []
        clean_list = clean_freq.to_a
        common_sum = clean_list.sum(){|a| a[1]}
        clean_list.each do |val|
          common_subseqs << val[0] if(val[1] > (common_sum * 0.05))

          #New and probably better:
          #common_subseqs << val[0] if(val[1] > (common_sum * 0.10) and val[1] > 3) #increased the thresholds as we were getting issues at some spots.
        end

        culm_time_i += Time.now - ts

        #realign to the options.
        final_options = [] #[alignment, orig, match]

        #possibly slightly slow
        ts = Time.now

        common_subseqs.each do |csubseq|
          #The $'s match super well together, preventing sequences from getting longer for better alignments.
          new_align = align_it(
            '$' + csubseq.gsub('-','') + '$',
            '$' + frag_seq.gsub('-','') + '$',
            0, 3)
          new_align[0].gsub!('$','')
          new_align[1].gsub!('$','')
          atmp = Alignment.new()
          atmp.std = Sequence.new(nuc: new_align[0])
          atmp.seq = Sequence.new(nuc: new_align[1])
          score = atmp.match_perc()
          final_options << [new_align, score] if(score > 0.70)
        end
        culm_time_b += Time.now - ts

        #choose best option
        best = final_options.max(){|a| a[1]}
        next if(!best) #No good alternate alignment?  Skip!

#        if(best)
#          log << "Indel:  #{indel},  range:  #{frag_rng},  clean_rng:  #{clean_rng}"
#          log << "Will replace:  #{frag_std}"
#          log << "               #{frag_seq}"
#          log << "With:          #{best[0][0]}"
#          log << "               #{best[0][1]}"
#          log << "Score:         #{best[1]}"
#        end


        #replace alignment if good, and also recalculate the insertion/deletion spots.
        #For the standard, DON'T replace, just move the gaps to the correct spots.
        std_frag = alignment.std.nucleotides[frag_rng].gsub('-','')

        #Moving the standard gaps to the correct spots.
        ts = Time.now
        fdex = 0
        best[0][0].each_char do |nuc|
          if(nuc == '-')
            if(fdex > std_frag.size())
              std_frag.insert(-1,'-')
            else
              std_frag.insert(fdex,'-')
            end
          end
          fdex += 1
        end
        culm_time_f += Time.now - ts


        needs_recalc = (alignment.std.nucleotides[frag_rng] != std_frag)
        ts = Time.now
        if(alignment.seq.nucleotides[frag_rng] != best[0][1] or needs_recalc)
          #log << "replacing..."
          #log << alignment.std.nucleotides[frag_rng.first - 3 .. frag_rng.last + 3]
          #log << alignment.seq.nucleotides[frag_rng.first - 3 .. frag_rng.last + 3]

          alignment.seq.nucleotides[frag_rng] = best[0][1]
          alignment.std.nucleotides[frag_rng] = std_frag
          alignment.cache_seq_clean = nil #important, since we modified the sequence.
          culm_time_c += Time.now - ts

          #log << "->"
          #log << alignment.std.nucleotides[frag_rng.first - 3 .. frag_rng.last + 3]
          #log << alignment.seq.nucleotides[frag_rng.first - 3 .. frag_rng.last + 3]
          #log << ""

          #Slowest!
          if(needs_recalc)
            ts = Time.now
            alignment.recalc_indels() #recalculate all the indels in the alignment.
            culm_time_d += Time.now - ts

            ts = Time.now
            alignment.seq_clean()
            culm_time_e += Time.now - ts
          end
        end

      end #end indels.each

      ts = Time.now
      alignment.cache_seq_clean =nil
      alignment.recalc_indels()
      culm_time_h += Time.now - ts
    end #end alignments.each
  end #end if(Settings['alignment-optimize'])

  #puts "Culm Time A: #{culm_time_a} seconds."
  #puts "Culm Time B: #{culm_time_b} seconds."
  #puts "Culm Time C: #{culm_time_c} seconds."
  #puts "Culm Time D: #{culm_time_d} seconds."
  #puts "Culm Time E: #{culm_time_e} seconds."
  #puts "Culm Time F: #{culm_time_f} seconds."
  #puts "Culm Time G: #{culm_time_g} seconds."
  #puts "Culm Time H: #{culm_time_h} seconds."
  #puts "Culm Time I: #{culm_time_i} seconds."

  #debug for conan
#  File.open("alignment-optimize.log", 'w') do |file|
#    log.each do |line|
#      file.puts line
#    end
#  end

end
