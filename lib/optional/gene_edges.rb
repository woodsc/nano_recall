require "#{File.dirname(__FILE__)}/../settings.rb"
require "#{File.dirname(__FILE__)}/../sample.rb"


#Seems like it actually makes things worse...
def trim_gene_edges(gene: )
  if(false and Settings['trim-gene-edges']) #disabled for now
    trim_start = 0
    trim_end = gene.hxb2_standard.nucleotides.size() - 1

    0.upto(gene.hxb2_standard.nucleotides.size() - 1) do |i|
      cnt = 0
      dash_cnt = 0
      gene.clades.each do |clade|
        clade.alignments.each do |alignment|
          cnt += 1
          dash_cnt += 1 if(alignment.seq_clean[i] == '-')
        end
      end
      if(dash_cnt.to_f / cnt.to_f > Settings['trim-edges-dash-threshold'])
        trim_start = i
      else
        break
      end
    end
    (gene.hxb2_standard.nucleotides.size() - 1).downto(0) do |i|
      cnt = 0
      dash_cnt = 0
      gene.clades.each do |clade|
        clade.alignments.each do |alignment|
          cnt += 1
          dash_cnt += 1 if(alignment.seq_clean[i] == '-')
        end
      end

      if(dash_cnt.to_f / cnt.to_f > Settings['trim-edges-dash-threshold'])
        trim_end = i - 1
      else
        break
      end
    end
    gene.trim_range = (trim_start .. trim_end)
  end #end if(Settings['trim-edges'])


end
