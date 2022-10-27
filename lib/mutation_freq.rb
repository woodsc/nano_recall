#Class that outputs mutation frequencies file.
require "#{File.dirname(__FILE__)}/settings.rb"
require "#{File.dirname(__FILE__)}/alignment.rb"
require "#{File.dirname(__FILE__)}/sequence.rb"
require "#{File.dirname(__FILE__)}/translations.rb"



class MutationFreq


  def MutationFreq.get_mutation_list(alignments: , hxb2_standard: , trim_range: )
    mut_list = [] # "Q12wt/T"
    freq = MutationFreq.get_aa_frequencies(alignments: alignments, trim_range: trim_range)

    freq.each_with_index() do |f_hsh, loc|
      muts = []
      hxb2_aa = nuc_to_aa(hxb2_standard.nucleotides[(loc * 3), 3])
      denominator = f_hsh.to_a.map(){|e| e[0] == '?' ? 0 : e[1] }.sum()

      f_hsh.to_a.each do |f|
        next if(f[0] == '?')
        perc = f[1].to_f / denominator.to_f
        if(perc > Settings['mix-threshold'])
          muts << (f[0] == hxb2_aa ? '!wt' : f[0] )
        end
      end

      if(muts != ['!wt'])
        mut_list << "#{hxb2_aa}#{loc + 1}" + muts.sort().join('/').gsub('!wt','wt')
      end
    end

    return mut_list
  end

  def MutationFreq.get_consensus_aa_array(alignments: , trim_range: )
    aa_consensus = []
    freq = MutationFreq.get_aa_frequencies(alignments: alignments, trim_range: trim_range)
    freq.each_with_index() do |f_hsh, loc|
      muts = []
      denominator = f_hsh.to_a.map(){|e| e[0] == '?' ? 0 : e[1] }.sum()

      f_hsh.to_a.each do |f|
        next if(f[0] == '?')
        perc = f[1].to_f / denominator.to_f
        if(perc > Settings['mix-threshold'])
          muts << f[0]
        end
      end

      aa_consensus << muts
    end

    return aa_consensus
  end

  def MutationFreq.get_consensus_aa_str(alignments: , trim_range: )
    #Not even sure how to do this?  Maybe just () the mixtures
  end

  def MutationFreq.get_consensus_nuc(alignments: , trim_range: )
    nuc_consensus = ''

    freqs = MutationFreq.get_nuc_frequencies(
      alignments: alignments, trim_range: trim_range, full_amino_only: true)


    freqs.each_with_index() do |f_hsh, loc|
      #puts f_hsh
      #exit
      nucs = []
      denominator = f_hsh.to_a.map(){|e| e[0] == '?' ? 0 : e[1] }.sum()

      n_max = [-1,'?']
      f_hsh.to_a.each do |f|
        #next if(f[0] == '?' or f[0] == '-' or f[0] == 'N')
        next if(f[0] == '?' or f[0] == 'N')
        perc = f[1].to_f / denominator.to_f
        if(perc > Settings['mix-threshold'])
          nucs << f[0] if(f[0] != '-')
          n_max = [perc, f[0]] if(perc > n_max[0])
        end
      end

      #puts "a: " + nucs.sort.inspect
      #puts "b: " + NUC_MIX.invert.inspect
      #puts "c: " + NUC_MIX.invert[nucs.sort()]
      if(n_max[1] == '-')  #We can't represent mixes of dels, so we chose del if its the greatest perc.
        nuc_consensus += '-'
      elsif(nucs != [])
        nuc_consensus += NUC_MIX.invert[nucs.sort()]
      else
        nuc_consensus += '-'
      end


    end

    return nuc_consensus
  end


  #Setting full_amino_only skips any nucleotides that don't form an amino acid.
  def MutationFreq.get_nuc_frequencies(alignments: , trim_range: , full_amino_only: false )
    freq = [] #array of {'A': 3, 'T': 2, 'C': 23, 'G':4}
    0.upto(alignments[0].seq_clean()[trim_range].size() - 1) do |i|
      freq[i] = {'A' => 0, 'T' => 0, 'C' => 0, 'G' => 0, '-' => 0, 'N' => 0, '?' => 0}
    end

    alignments.each do |alignment|
      seq_clean = alignment.seq_clean()[trim_range]

      #we need to only include the start to the end of the sequence, trimming edges.
      n_st = 0
      n_end = seq_clean.size() - 1
      n_st += 1 while(n_st < seq_clean.size() and seq_clean[n_st] == '-')
      n_end -= 1 while(n_end >= 1 and seq_clean[n_end] == '-')

      if(full_amino_only)
        0.upto((seq_clean.size() / 3).to_i - 1) do |i|
          if(!seq_clean[i*3, 3].include?('-') or (i >= n_st && i <= n_end && seq_clean[i*3, 3] == '---'))
            freq[i*3 + 0][seq_clean[i*3 + 0]] += 1
            freq[i*3 + 1][seq_clean[i*3 + 1]] += 1
            freq[i*3 + 2][seq_clean[i*3 + 2]] += 1
          end
        end
      else
        0.upto(seq_clean.size() - 1) do |i|
          if(i >= n_st && i <= n_end) #prevents sequence edge dashes from being included
            freq[i][seq_clean[i]] += 1
          end
        end
      end
    end
    return freq
  end

  #returns an array like:  NUC_LOC,AA_LOC,CNT,PERC,SIZE,NUC,AA,VALID
  def MutationFreq.get_ins_frequencies(alignments: )
    list = []

    alignments.each do |alignment|
      alignment.insertions.each do |ins|
        ins_loc = alignment.loc_ref(ins.begin)

        ins_exists = list.find(){|a| a[0] == ins_loc and a[5] == alignment.seq.nucleotides[ins]}
        if(ins_exists)
          ins_exists[2] += 1
        else
          list.push([
            ins_loc,
            (ins_loc / 3),
            1, #cnt
            nil,
            ins.size(),
            alignment.seq.nucleotides[ins],
            nuc_to_aa_string(alignment.seq.nucleotides[ins]),
            (ins.size() % 3 == 0 and ins_loc % 3 == 0 )  #validity (size+framealign)
          ])
        end
      end
    end

    #add the percent
    list.each do |row|
      row[3] = ((row[2].to_f / alignments.size().to_f) * 100).round(1)
    end

    #now we order by loc, count
    list.sort!() do |a,b|
      if(a[0] == b[0])
        b[2] <=> a[2]
      else
        a[0] <=> b[0]
      end
    end

    return list
  end

  def MutationFreq.get_aa_frequencies(alignments: , trim_range: )
    return [] if(alignments.size() == 0)

    freq = []

    0.upto((alignments[0].seq_clean()[trim_range].size() / 3) - 1) do |i|
      freq[i] = {}
      freq[i].default = 0
    end

    alignments.each do |alignment|
      seq = alignment.seq_clean()[trim_range]
      aa_seq = ''
      0.upto((seq.size() / 3).to_i - 1) do |i|
        aa = nuc_to_aa(seq[i * 3, 3])
        aa_seq += (aa.nil? ? '?' : aa)
        freq[i][(aa.nil? ? '?' : aa)] += 1 if(freq[i])
      end

    end
    return freq
  end

  def MutationFreq.save_nuc(filename: , alignments: , hxb2_standard: , trim_range: )
    freq = MutationFreq.get_nuc_frequencies(alignments: alignments, trim_range: trim_range)

    #inserts
    ins_freq = {}
    ins_freq.default = 0
    alignments.each do |alignment|
      alignment.insertions.each do |ins|
        ins_freq["#{alignment.loc_ref(ins.first)}-#{alignment.loc_ref(ins.first) + ins.size()}:#{alignment.seq.nucleotides[ins]}"] += 1
      end
    end

    File.open(filename, 'w') do |file|
      row = ['Loc','HXB2','CALL','SECONDARY']
      "ACGT-N".split('').each do |nuc|
        row << "#{nuc} count"
      end
      row << ''
      "ACGT-N".split('').each do |nuc|
        row << "#{nuc} perc"
      end

      file.puts row.join(',')

      freq.each_with_index do |fr, dex|
        freq_order = fr.to_a.sort(){|a,b| b[1] <=> a[1]}
        secondary_perc = freq_order[1][1].to_f / freq_order.map(){|a| a[1]}.sum().to_f
        row = [dex+1,
          hxb2_standard.nucleotides[dex],
          freq_order[0][0],
          (secondary_perc > Settings['mix-threshold'] ? freq_order[1][0] : '')
        ]
        "ACGT-N".split('').each do |nuc|
          row << fr[nuc]
        end
        row << ''
        "ACGT-N".split('').each do |nuc|
          row << ((fr[nuc].to_f / alignments.size().to_f) * 100).round(2)
        end
        file.puts row.join(',')
      end

      file.puts
      file.puts "Insertion frequencies:"
      file.puts "Loc,Size,HXB2,Insert,Count,Percent"
      ins_freq.to_a.sort(){|a,b| a[0].split('-')[0].to_i <=> b[0].split('-')[0].to_i}.each do |val|
        fr_key, fr_value = val[0], val[1]
        loc = fr_key.split('-')[0].to_i
        size = fr_key.split('-')[1].split(':')[0].to_i - loc
        insert = fr_key.split(':')[1]
        file.puts "#{loc+1},#{size},#{hxb2_standard.nucleotides[loc]},#{insert},#{fr_value},#{((fr_value.to_f / alignments.size().to_f) * 100).round(2)}"
      end
    end
  end

  def MutationFreq.save_aa(filename: , alignments: , hxb2_standard: , trim_range: )
    aa_columns = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','*','-','ins']
    freq = MutationFreq.get_aa_frequencies(alignments: alignments, trim_range: trim_range)

    File.open(filename, 'w') do |file|
      row = ['LOC','HXB2', 'CALL','SECONDARY', 'PERCENTS'] + aa_columns.map(){|a| (a == '-' ? 'del' : a ) + " cnt"}

      file.puts row.join(',')

      freq.each_with_index do |f_hsh, loc|
        row = [loc + 1, nuc_to_aa(hxb2_standard.nucleotides[loc * 3, 3]), nil, nil, nil ]
        sorted_percs = f_hsh.to_a.sort(){|a,b| b[1] <=> a[1] }.select(){|e| e[0] != '?'}
        denominator = sorted_percs.map(){|e| e[0] == '?' ? 0 : e[1]  }.sum().to_f
        if(denominator != 0)
          row[2] = sorted_percs[0][0] if(sorted_percs[0])
          row[3] = sorted_percs.find_all() {|e|
            (e[0] != row[2] and e[1].to_f / denominator > Settings['mix-threshold'])
          }.map(){|e| e[0] }.sort().join('')

          row[4] = sorted_percs.
            find_all(){|e| ((e[1].to_f / denominator) * 100).to_i > 0}.
            map(){|e| e[0] + ((e[1].to_f / denominator) * 100).to_i.to_s + '%'  }.
            join(' ')
        end

        aa_columns.each_with_index do |aa_col, aa_dex|
          row[5 + aa_dex] = f_hsh[aa_col] ? f_hsh[aa_col] : 0
        end

        ins_col = 5 + aa_columns.index('ins')
        row[ins_col] = 0
        alignments.each do |alignment|
          alignment.insertions.each do |ins|
            ins_loc = alignment.loc_ref(ins.begin)
            if((ins_loc / 3) == loc + 1 and ins.size >= 3 and ins.size % 3 == 0)
              row[ins_col] += 1
              #puts "#{alignment.seq.id}"
              #puts "#{ins_loc} #{ins_loc / 3} #{alignment.seq.nucleotides[ins]} -> #{nuc_to_aa_string(alignment.seq.nucleotides[ins])}"
              break
            end
          end
        end


        file.puts row.join(',')
      end

    end

  end

end
