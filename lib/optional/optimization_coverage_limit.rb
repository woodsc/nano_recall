require "#{File.dirname(__FILE__)}/../settings.rb"
require "#{File.dirname(__FILE__)}/../sample.rb"


def optimization_coverage_limit(gene: , alignment: , region_def:)
  if(Settings['optimization-coverage-limit'])
    #update amino acid coverage.
    seq = alignment.seq_clean()
    aa_seq = ''
    0.upto((seq.size() / 3).to_i - 1) do |i|
      aa = nuc_to_aa(seq[i * 3, 3])
      aa_seq += (aa.nil? ? '?' : aa)
    end

    #Clear out the edges, as its like atg,a--,---,---...
    aa_seq.gsub!(/^(\-+\?+)/){|match| '-' * match.size()}
    aa_seq.gsub!(/(\?+\-+)$/){|match| '-' * match.size()}

    0.upto(aa_seq.size() - 1) do |i|
      gene.aa_coverage[i] += 1 if(gene.aa_coverage[i] != nil and aa_seq[i] != nil and aa_seq[i] != '?')
    end


    all_coverage = true
    0.upto(gene.aa_coverage.size() - 1) do |i|
      next if(i < region_def[0] or i > region_def[1])

      if(gene.aa_coverage[i] < Settings['optimization-coverage-limit-target'])
        all_coverage = false
        break
      end
    end


    if(all_coverage)
    #if(gene.aa_coverage.all?(){|a| a >= Settings['optimization-coverage-limit-target'] } )
      puts "#{gene.name} - Target coverage reached (min cov: #{gene.region_aa_coverage.min}  max cov: #{gene.region_aa_coverage.max()})."
      return true #quit_gene_early = true
    end


  end #end if(Settings['optimization-coverage-limit'])
  return false
end
