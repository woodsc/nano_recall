require "#{File.dirname(__FILE__)}/settings.rb"
require "#{File.dirname(__FILE__)}/asi_algorithm.rb"

class ResistanceReport
  attr_accessor :label, :pr_result, :rt_result, :int_result
  attr_accessor :pr_aa, :rt_aa, :int_aa
  attr_accessor :pr_aa_hxb2, :rt_aa_hxb2, :int_aa_hxb2
  attr_accessor :algorithm

  def initialize(label: '', algorithm: nil)
    @label = label
    @algorithm = algorithm.nil? ? AsiAlgorithm.new("config/HIVDB_9.0.xml") : algorithm
    @pr_result = nil
    @rt_result = nil
    @int_result = nil
  end


  def add_aa(aa_seq: , gene: , aa_hxb2: )
    if(gene == 'pr')
      @pr_aa = trim_aa(aa_seq)
      @pr_aa_hxb2 = aa_hxb2
      @pr_result = @algorithm.interpret(@pr_aa, 'PR')
    elsif(gene == 'rt')
      @rt_aa = trim_aa(aa_seq)
      @rt_aa_hxb2 = aa_hxb2
      @rt_result = @algorithm.interpret(@rt_aa, 'RT')
    elsif(gene == 'int')
      @int_aa = trim_aa(aa_seq)
      @int_aa_hxb2 = aa_hxb2
      @int_result = @algorithm.interpret(@int_aa, 'IN')
    end

  end

  #removes dashes from begining and end (replaces with empty)
  def trim_aa(aa_seq)
    trimmed = aa_seq.clone()
    0.upto(trimmed.size() - 1) do |i|
      if(trimmed[i] == ['-'])
        trimmed[i] = []
      else
        break
      end
    end
    (trimmed.size() - 1).downto(0) do |i|
      if(trimmed[i] == ['-'])
        trimmed[i] = []
      else
        break
      end
    end
    return trimmed
  end

  #creates a mutation list based on hxb2
  def mut_list(aa, hxb2)
    list = []
    0.upto(aa.size() - 1) do |i|
      if(aa[i].size() == 1 and aa[i].first == hxb2[i])
        #pass
      elsif(aa[i] == [])
        #pass
      else
        tmp = aa[i].sort()
        if(aa[i].include?(hxb2[i]))
          tmp = ['wt'] + (aa[i] - [hxb2[i]]).sort()
        end
        list << "#{hxb2[i]}#{i+1}#{tmp.join('')}"
      end
    end
    return list
  end

  def key_mut_list(region, aa, hxb2)
    list = []
    mut_list = mut_list(aa, hxb2)
    result = nil
    if(region == 'pr')
      result = @pr_result
    elsif(region == 'rt')
      result = @rt_result
    elsif(region == 'int')
      result = @int_result
    end

    mut_list.each do |mut|
      mut =~ /^\D(\d+)\D+$/
      loc = $1.to_i
      tmp = result.mutations.find() {|a| a[0] == loc }
      list << mut if(tmp)
    end

    #result.mutations.each do |mut|
    #  tmp = mut_list.find(){|a| a =~ /^\D(\d+)\D+$/;$1.to_i == mut[0] }
    #  list << tmp
    #end

    return list
  end


  def text_report()
    text = "Resistance report - #{@label}\n"
    text += "For research purposes only\n\n"

    #longest drug name, used for padding strings.
    pad = ((@pr_result ?  @pr_result.drugs : []) +
      (@rt_result ?  @rt_result.drugs : []) +
      (@int_result ?  @int_result.drugs : [])).
      map(){|a| a.name.size()}.max()
    pad = (pad.nil? ? 2 : pad + 2)

    if(@pr_result)
      text += "Protease Inhibitors:\n"
      @pr_result.drugs.each do |drug|
        call = @algorithm.level_def.find(){|e| drug.level.to_s == e[0]}
        text += "  #{"#{drug.name.capitalize()}:".ljust(pad, ' ')}  #{call[1]}\n"
      end
      text += "Key Protease Mutations:\n"
      text += "  #{key_mut_list('pr', @pr_aa, @pr_aa_hxb2).join(' ')}\n\n"
      text += "All Protease Mutations:\n"
      text += "  #{mut_list(@pr_aa, @pr_aa_hxb2).join(' ')}\n\n"
    end

    if(@rt_result)
      text += "Nucleoside Reverse Transcriptase Inhibitors:\n"
      @rt_result.drugs.find_all(){|d| d.drug_class == 'NRTI'}.each do |drug|
        call = @algorithm.level_def.find(){|e| drug.level.to_s == e[0]}
        text += "  #{"#{drug.name.capitalize()}:".ljust(pad, ' ')}  #{call[1]}\n"
      end

      text += "Non-nucleoside Reverse Transcriptase Inhibitors:\n"
      @rt_result.drugs.find_all(){|d| d.drug_class == 'NNRTI'}.each do |drug|
        call = @algorithm.level_def.find(){|e| drug.level.to_s == e[0]}
        text += "  #{"#{drug.name.capitalize()}:".ljust(pad, ' ')}  #{call[1]}\n"
      end
      text += "Key Reverse Transcriptase Mutations:\n"
      text += "  #{key_mut_list('rt', @rt_aa, @rt_aa_hxb2).join(' ')}\n\n"
      text += "All Reverse Transcriptase Mutations:\n"
      text += "  #{mut_list(@rt_aa, @rt_aa_hxb2).join(' ')}\n\n"
    end

    if(@int_result)
      text += "Integrase Strand Transfer Inhibitors:\n"
      @int_result.drugs.each do |drug|
        call = @algorithm.level_def.find(){|e| drug.level.to_s == e[0]}
        text += "  #{"#{drug.name.capitalize()}:".ljust(pad, ' ')}  #{call[1]}\n"
      end
      text += "Key Integrase Mutations:\n"
      text += "  #{key_mut_list('int', @int_aa, @int_aa_hxb2).join(' ')}\n\n"
      text += "All Integrase Mutations:\n"
      text += "  #{mut_list(@int_aa, @int_aa_hxb2).join(' ')}\n\n"
    end

    if(@pr_result)
      text += "Protease Mutation Comments:\n"
      @pr_result.mutation_comments.each do |comment|
        text += "  #{comment}\n"
      end
    end

    if(@rt_result)
      text += "Reverse Transcriptase Mutation Comments:\n"
      @rt_result.mutation_comments.each do |comment|
        text += "  #{comment}\n"
      end
    end

    if(@int_result)
      text += "Integrase Mutation Comments:\n"
      @int_result.mutation_comments.each do |comment|
        text += "  #{comment}\n"
      end
    end

    mix_perc = (Settings['mix-threshold'] * 100).round(1)
    text += "\nMutations were included if they made up #{mix_perc}% of the sequences analyzed.\n"

    return text
  end


end
