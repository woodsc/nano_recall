=begin
Sweeps close by inserts into a single big insertion, and then tries to frame
align the insertion to the best spot.  If the sequence has a LOT of single
base insertions, it may have poor results...
=end

def sweep_inserts(alignment: )
  #puts "#{alignment.seq.id}  Original score:  #{alignment.match_perc()}"
  #puts alignment.std.nucleotides[alignment.trim_start,60]
  #puts alignment.seq.nucleotides[alignment.trim_start,60]
  #puts alignment.insertions.inspect

  #step 1, identify insertions to sweep together
  state = :start
  clusters = []
  combo = []
  ins_dex = 0

  while(ins_dex < alignment.insertions.size() + 1) #we go one further for cleanup purposes.
    #keep chaining nearby insertions
    cur_ins = alignment.insertions[ins_dex]

    #fix for end of insert list chores.
    if(cur_ins == nil && state == :chain) #end of inserts.
      state = :clustered
    end

    if(state == :start)
      combo << cur_ins
      ins_dex += 1
      state = :chain
    elsif(state == :chain)
      if(cur_ins.begin - combo.last.end <= 6)
        combo << cur_ins
        ins_dex += 1
      else
        state = :clustered
      end
    elsif(state == :clustered)
      #create a list of insertions to cluster.
      if(combo.map(){|a| a.size()}.sum() >= 3)
        clusters << combo
        combo = []
        state = :start
        ins_dex += 1
      else
        combo = []
        state = :start
      end
    end
  end


  #now, we try to put the insertions into frame-aligned, high scoring locations.
  clusters.each do |cluster|
    options = []
    region = ([cluster.first.begin - 6, 0].max() .. [cluster.last.end + 6, alignment.std.nucleotides.size() - 1].min())
    ins_size = cluster.map(){|a| a.size()}.sum()

    #create list of candidate alignments.
    region.each do |dex|
      next if(alignment.std.nucleotides[dex] == '-') #avoids a possible bug in loc_ref???
      if(alignment.loc_ref(dex) % 3 == 0 and (dex <= region.end - ins_size)) #frame aligned.
        std = alignment.std.nucleotides[region]

        cluster.reverse().each do |ins|
          #using O as a placeholder, so we can delete it later without messing up sequence size
          std[(ins.begin - region.begin) .. (ins.end - region.begin)] = 'O' * ins.size()
        end
        std.insert(dex - region.begin, '-' * ins_size)
        std.gsub!('O','')
        options << [std]
      end
    end

    next if(options.empty?())
    
    #choose frame aligned area with highest match score.
    options.each do |opt|
      tmp = Alignment.new
      tmp.std = Sequence.new(nuc: opt[0])
      tmp.seq = Sequence.new(nuc: alignment.seq.nucleotides[region], id: 'tmp')
      opt.push(tmp.match_perc)
    end

    pick = options.max(){|a,b| a[1] <=> b[1]}[0]
    alignment.std.nucleotides[region] = pick
    alignment.recalc_indels()
  end

end
