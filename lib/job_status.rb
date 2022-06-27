require "#{File.dirname(__FILE__)}/settings.rb"
require "#{File.dirname(__FILE__)}/sample.rb"
require "#{File.dirname(__FILE__)}/resistance_report"
require 'json'


def job_setup(filename, genes)
  job = {
    file: filename,
    progress: 0, #0-100 percent complete, use the min of the cutoff/num_seqs to calc
    status: 'ip',
    num_sequences: nil,
    resistance_report: nil,
    date_complete: nil,
    start_ts: Time.now,
    end_ts: nil,
    update_ts: Time.now - 10,
    genes: {},  #which contains name, num_analyzed, mutations, subtype percs.
  }

  genes.each do |gene_name|
    job[:genes][gene_name] = {
      name: gene_name,
      num_analyzed: 0,
      mutations: [],
      subtypes: [],
    }
  end
  return job
end

def job_add_num_sequences(job, num)
  job[:num_sequences] = num
end

def job_update(dex, job, sample, force=false)
  update_interval = Settings['jobfile-update-interval']
  #if 15 seconds have gone by, update.
  if((job[:update_ts] + update_interval) < Time.now() or force)
    cov_limit = Settings['optimization-coverage-limit-target']
    completion = {}
    if(Settings['optimization-coverage-limit'] and cov_limit.to_i < job[:num_sequences])
      sample.genes.each do |sample_gene|
        cov = sample_gene.aa_coverage.min()
        completion[sample_gene.name] = [cov, cov_limit.to_i]
      end
    else
      sample.genes.each do |sample_gene|
        completion[sample_gene.name] = [dex, job[:num_sequences]]
      end
    end
    #puts "Analyzed sequence #{dex} of #{job[:num_sequences]}"
    min_perc = 100
    sample.genes.each do |sample_gene|
      g_comp = completion[sample_gene.name]
      perc = ((g_comp[0].to_f / g_comp[1].to_f) * 100).to_i
      min_perc = perc if(perc < min_perc)
      #puts "#{sample_gene.name}: #{g_comp[0]} of #{g_comp[1]}   (#{perc}%)."
    end
    job[:progress] = min_perc
    #puts "Total progress: #{min_perc}%"

    #in-progress report
    report = ResistanceReport.new(label: sample.id)
    sample.genes.each do |gene|
      gene_alignments = []
      gene.clades.each {|clade| gene_alignments += clade.alignments}

      next if(gene_alignments.size < 10) #not enough data for even prelim report

      #Add data to resistance report
      aa_consensus = MutationFreq.get_consensus_aa_array(
        alignments: gene_alignments, trim_range: gene.trim_range)

      report.add_aa(aa_seq: aa_consensus, gene: gene.name,
        aa_hxb2: nuc_to_aa_array(gene.hxb2_standard.nucleotides) )

      job[:genes][gene.name][:num_analyzed] = gene_alignments.size()
      job[:genes][gene.name][:subtypes] = []
      gene_seq_cnt = gene.clades.map(){|c| c.alignments.size()}.sum()
      gene.clades.sort(){|a,b| b.alignments.size() <=> a.alignments.size() }.each do |clade|
        job[:genes][gene.name][:subtypes] << [
          clade.id.split(' ')[1],
          clade.alignments.size(),
          ((clade.alignments.size().to_f / gene_seq_cnt.to_f) * 100).to_i
        ]
      end
    end
    job_update_report(job, report)

    job[:update_ts] = Time.now()
    #pp job[:resistance_report]
    #puts
    return true
  end
  return false
end

#updates the job based on the report object
def job_update_report(job, report)
  #in-progress report
  if(report.pr_aa)
    job[:genes]['pr'][:mutations] = report.mut_list(report.pr_aa, report.pr_aa_hxb2)
  end
  if(report.rt_aa)
    job[:genes]['rt'][:mutations] = report.mut_list(report.rt_aa, report.rt_aa_hxb2)
  end
  if(report.int_aa)
    job[:genes]['int'][:mutations] = report.mut_list(report.int_aa, report.int_aa_hxb2)
  end

  rep = job[:resistance_report] = {}
  rep[:alg_version] = "#{report.algorithm.alg_name} #{report.algorithm.alg_version}"
  rep[:pr] = {}
  rep[:rt] = {}
  rep[:int] = {}
  rep[:pr][:drugs] = []
  rep[:rt][:drugs] = []
  rep[:int][:drugs] = []
  rep[:pr][:key_mutations] = []
  rep[:rt][:key_mutations] = []
  rep[:int][:key_mutations] = []
  rep[:pr][:comments] = []
  rep[:rt][:comments] = []
  rep[:int][:comments] = []

  if(report.pr_result)
    report.pr_result.drugs.each do |drug|
      call = report.algorithm.level_def.find(){|e| drug.level.to_s == e[0]}
      rep[:pr][:drugs] << {
        drug: drug.name,
        code: drug.code,
        drug_class: drug.drug_class,
        drug_level: drug.level,
        call: call[1],
      }
    end
    rep[:pr][:key_mutations] = report.key_mut_list('pr', report.pr_aa, report.pr_aa_hxb2)
    rep[:pr][:comments] = report.pr_result.mutation_comments
  end
  if(report.rt_result)
    report.rt_result.drugs.each do |drug|
      call = report.algorithm.level_def.find(){|e| drug.level.to_s == e[0]}
      rep[:rt][:drugs] << {
        drug: drug.name,
        code: drug.code,
        drug_class: drug.drug_class,
        drug_level: drug.level,
        call: call[1],
      }
    end
    rep[:rt][:key_mutations] = report.key_mut_list('rt', report.rt_aa, report.rt_aa_hxb2)
    rep[:rt][:comments] = report.rt_result.mutation_comments
  end
  if(report.int_result)
    report.int_result.drugs.each do |drug|
      call = report.algorithm.level_def.find(){|e| drug.level.to_s == e[0]}
      rep[:int][:drugs] << {
        drug: drug.name,
        code: drug.code,
        drug_class: drug.drug_class,
        drug_level: drug.level,
        call: call[1],
      }
    end
    rep[:int][:key_mutations] = report.key_mut_list('int', report.int_aa, report.int_aa_hxb2)
    rep[:int][:comments] = report.int_result.mutation_comments
  end
end

def job_complete(job, report)
  job[:status] = 'success'
  job[:progress] = 100
  job[:end_ts] = Time.now

  #update report
  job_update_report(job, report)
end

def job_error(job, error)
  job[:status] = 'fail'
  job[:progress] = 100
  job[:end_ts] = Time.now
  job[:error] = error
end

def job_save(job)
  if(job[:file])
    File.open(job[:file], 'w') do |file|
      file.puts job.to_json()
    end
  end
end


def job_zip(job, path)
  begin
    FileUtils.rm("#{path}.zip")
  rescue
  end
  
  if(job[:file])
    system("zip -j #{path}.zip #{path}\.*  #{job[:file]}")
  else
    system("zip -j #{path}.zip #{path}\.*")
  end
end
