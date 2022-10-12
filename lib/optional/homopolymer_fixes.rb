require "#{File.dirname(__FILE__)}/../settings.rb"
require "#{File.dirname(__FILE__)}/../sample.rb"

=begin
1)  Searches the standard for homopolymers (3 or greater)
For each homopolymer,
2)  Realign to the most common sequence in this region.  Replace if better match.
3)  Fill in gaps in simple homopolymer regions
=end

#I think this is right, but we need to be careful and test and clean.
def homopolymer_fixes_get_homopolymers(gene: , apad: )
  #first pass
  homopolymers = []
  st = -1
  hnuc = ''
  gene.hxb2_standard.nucleotides.split('').each_with_index do |nuc, i|
    if(st == -1)
      st = i
      hnuc = nuc
    elsif(hnuc != nuc) #end homopolyer
      if(i - st >= 3)
        homopolymers << [st, (i - st), gene.hxb2_standard.nucleotides[st, (i - st)]]
      end
      st = i
      hnuc = nuc
    end
  end

  #second pass, find most common sequences and correct homopolymer loc/size
  #Add in most common sequence.
  homopolymers.each do |hp|
    clean_list = []
    clean_freq = []

    gene.clades.each do |clade|
      clade.alignments.each do |alignment|
        subseq = alignment.seq_clean[hp[0] - apad, hp[1] + apad*2]
        clean_list << subseq if(!subseq.include?('-'))
      end
    end

    clean_list.uniq.each do |seq|
      clean_freq << [seq, clean_list.find_all(){|e| e == seq}.size()]
    end
    common_subseq = clean_freq.sort(){|a,b| b[1] <=> a[1]}.first()[0]
    #correct homopolymer size/region
    kill = true
    xst = -1
    xhnuc = ''
    common_subseq.split('').each_with_index do |nuc, i|
      if(xst == -1)
        xst = i
        xhnuc = nuc
      elsif(xhnuc != nuc) #end homopolyer
        if(i - xst >= 3)
          if((xst - apad).abs() <= 2) #make sure we are at most 2 away from original loc.
            hp_orig = hp.clone()
            hp[0] += (xst - apad)
            hp[1] = (i - xst)
            hp[2] = common_subseq[xst, (i - xst)]
            kill = false
            break
          end
        end
        xst = i
        xhnuc = nuc
      end
    end
    hp[0] = nil if(kill)
  end

  homopolymers.delete_if(){|a| a[0] == nil}
  return homopolymers
end

#NEed to simplify and extract bits out to methods.
def homopolymer_fixes(gene: )
  if(Settings['homopolymer-fixes'])
    homopolymers = []
    apad = 6 #alignment_pad

    #Correct homopolymers based on the most common sequence,
    homopolymers = homopolymer_fixes_get_homopolymers(gene: gene , apad: apad)
    #pp homopolymers

    homopolymers.each do |hp|
      clean_list = []
      clean_freq = []

      substd = gene.hxb2_standard.nucleotides[hp[0] - apad, hp[1] + apad*2]


      #get the most common subsequence
      gene.clades.each do |clade|
        clade.alignments.each do |alignment|
          subseq = alignment.seq_clean[hp[0] - apad, hp[1] + apad*2]
          clean_list << subseq if(!subseq.include?('-'))
        end
      end
      clean_list.uniq.each do |seq|
        clean_freq << [seq, clean_list.find_all(){|e| e == seq}.size()]
      end
      common_subseq = clean_freq.sort(){|a,b| b[1] <=> a[1]}.first()[0]
      #puts
      #puts substd + '  (standard  ) '
      #puts common_subseq + '  (common seq) ' + hp.to_s

      #main work
      cum_score_orig = 0
      cum_score_new = 0
      gene.clades.each do |clade|
        clade.alignments.each do |alignment|
          subseq = alignment.seq_clean[hp[0] - apad, hp[1] + apad*2]
          new_align = align_it(common_subseq, subseq.gsub('-',''), 3, 1)
          if(new_align[0].size() > subseq.size()) #No good, skip
            new_align = [common_subseq, subseq]
            #puts 'bad'
          end

          matches_orig = 0
          matches_new = 0
          0.upto(subseq.size() - 1) do |i|
            matches_orig += 1 if(subseq[i] == common_subseq[i])
            matches_new += 1 if(new_align[1][i] == common_subseq[i])
          end
          cum_score_orig += matches_orig
          cum_score_new += matches_new


          if(matches_new > matches_orig and new_align[1].size() == subseq.size())
            #puts "#{subseq}  ->  #{new_align[1]}   changed:   #{matches_orig} #{matches_new}"
            #update alignment
            alignment.seq_clean[hp[0] - apad, hp[1] + apad*2] = new_align[1]
          else
            #puts "#{subseq}  ->  #{new_align[1]}   unchanged: #{matches_orig} #{matches_new}"
          end

          #Optionally fill in homopolymer gaps
          if(Settings['homopolymer-fixes-fill-gaps'])
            filled = false
            asc = "" + alignment.seq_clean
            hp[0].upto(hp[0] + hp[1] - 1) do |i|
              if(alignment.seq_clean[i] == '-' and
                common_subseq[i + apad - hp[0]] == hp[2][0] and
                alignment.seq_clean[hp[0]-1] != '-' and
                alignment.seq_clean[hp[0]+hp[1]] != '-')

                alignment.seq_clean[i] = hp[2][0]
                filled = true
              end
            end
            #puts "#{asc[hp[0] - apad, hp[1] + apad*2]}  ->  #{alignment.seq_clean[hp[0] - apad, hp[1] + apad*2]}   (Filled)" if(filled)
          end

        end
      end
      #puts "cum score change #{cum_score_orig} -> #{cum_score_new}"


    end

  end
end



#new

#experimental: shortens suspicious dual-homopolymers.
def reduce_dual_homopolymers(alignment: )
  if(Settings['reduce-dual-homopolymers'])
    #Attempts to shorten suspicious dual-homopolymers, EG  AAAAGGGG
    homopolymers = []
    st = -1
    hnuc = ''

    #First, find all the homopolymer
    alignment.trim_start.upto(alignment.trim_end) do |i|
      nuc = alignment.seq.nucleotides[i]
      if(st == -1)
        st = i
        hnuc = nuc
      elsif(hnuc != nuc) #end homopolyer
        if(i - st >= 3 and hnuc != '-')
          homopolymers << [st, (i - st), alignment.seq.nucleotides[st, (i - st)]]
        end
        st = i
        hnuc = nuc
      end
    end

    #Find two in a row, where the standard sequence is altered one way or another.  EG:
    #AAAGGGG
    #AAAAGGG
    #Would end up censoring as AAA-GGG
    #Its a bit complex though, so be careful.
    prev_hp = nil
    homopolymers.each do |hp|
      if(prev_hp == nil)
        prev_hp = hp
      elsif(prev_hp[0] + prev_hp[1] == hp[0] and hp[01] + prev_hp[1] >= 7)
        if(alignment.std.nucleotides[hp[0]] == prev_hp[2][0] and hp[1] > 3)
          alignment.seq.nucleotides[hp[0]] = '-'
          #tmp_seq = alignment.seq.nucleotides[prev_hp[0]-3, prev_hp[1] + hp[1]+6]
          #puts "A:  Should be:  #{tmp_seq}" if(final_alignment.seq.id == '3b8785b1-5816-4bbb-8b45-9d4a89dd4f39')
        elsif(alignment.std.nucleotides[hp[0] - 1] == hp[2][0] and prev_hp[1] > 3)
          alignment.seq.nucleotides[hp[0] - 1] = '-'
          #tmp_seq = alignment.seq.nucleotides[prev_hp[0]-3, prev_hp[1] + hp[1]+6]
          #puts "B:  Should be:  #{tmp_seq}" if(final_alignment.seq.id == '3b8785b1-5816-4bbb-8b45-9d4a89dd4f39')
        end
      end
      prev_hp = hp
    end

  end
end
