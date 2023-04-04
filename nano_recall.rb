require 'optparse'
require 'pp'
require 'fileutils'


require './lib/settings'
require './lib/job_status'
require './lib/fasta'
require './lib/fastq'
require './lib/sequence'
require './lib/alignment'
require './lib/sample'
require './lib/mutation_freq'
require './lib/resistance_report'
require './lib/alignment/_alignment' #might get rid of this and use it in alignment.rb

require './lib/optional/gene_edges'
require './lib/optional/homopolymer_fixes'
require './lib/optional/alignment_optimize'
require './lib/optional/optimization_coverage_limit'
require './lib/optional/reject_poor_aligned_ends'
require './lib/optional/sweep_inserts'


options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: nano_recall.rb [options]"

  opts.on("-h", "--help", "Prints help") do |v|
    options[:help] = v
    puts opts
  end
#  opts.on("-c", "--command COMMAND", "Runs a command [call,zip]") do |v|
#    options[:command] = v
    #check that it is valid
#  end
  opts.on("-i", "--input FILE", "Input file") do |v|
    options[:input] = v
    #check that it is a valid file
  end
  opts.on("--batch-input", "--batch-input FOLDER", "Input folder") do |v|
    options[:batch_input] = v
    #check that it is a valid file
  end
  opts.on("-o", "--output FILE", "Output file") do |v|
    options[:output] = v
    #check that it is valid output file?
  end
  opts.on("-z", "--zip", "Zip output files") do |v|
    options[:zip] = v
    #check that it is valid output file?
  end
  opts.on("--batch-output", "--batch-output FOLDER", "Output folder") do |v|
    options[:batch_output] = v
    #check that it is valid output file?
  end
  opts.on("--use-sampleids", "--use-sampleids", "Use sample ids for batch output") do |v|
    #enable sample name output
    options[:use_sampleids] = true
  end
#  opts.on("-p", "--project PROJECT", "Project to use") do |v|
#    options[:project] = v
#  end
  opts.on("-d", "--inspect ID", "Inspect a specific ID (DEBUGING)") do |v|
    options[:inspect] = v
  end
  opts.on("-j", "--job FILE", "Choose file to hold ip results.  (batch processing not supported)") do |v|
    options[:jobfile] = v
    #check that it is valid output file?
  end
end.parse!

batch_files = []
stop = false

if(options[:help])
  stop = true
else
  #if(!options[:command])
  #  puts "Missing required option --command"
  #  stop = true
  #end
  if(options[:batch_output] or options[:batch_input])
    #process everything in the batch_input folder.
    options[:batch_input].gsub!(/[\/\\]$/,'')
    files = Dir["#{options[:batch_input]}/*.gz"] + Dir["#{options[:batch_input]}/*.fastq"]
    files.each do |file|
      out = "#{options[:batch_output]}/#{file.split('/').last.gsub('.gz','').gsub('.fastq','')}"
      batch_files << {input: file, output: out}
    end
  else
    if(!options[:input])
      puts "Missing required option --input"
      stop = true
    end
    if(!options[:output])
      puts "Missing required option --output"
      stop = true
    end
    batch_files << {input: options[:input], output: options[:output]}
  end

end
if(stop)
  exit(1)
end



Settings.load("./config/settings.txt")
quit_early = false #if true, we stop processing early.
if(Settings['debug'])
  pp options
end


#Load projects.
projects = nil
standards = load_fasta(filename: "config/reference_standards.fas", reference: true)

gene_list = ['pr', 'rt', 'int']
gene_list = Settings['debug-genes'].split(',')

#Create Folders if needed.
if(options[:jobfile] =~ /^(.+)\/[^\/]+$/)
  FileUtils.mkdir_p($1)
end
if(options[:output] =~ /^(.+)\/[^\/]+$/)
  FileUtils.mkdir_p($1)
end

#Set up job
job = job_setup(options[:jobfile], gene_list)

stats = {
  no_match_cnt: 0,
  reverse_cnt: 0,
  forward_cnt: 0,
  total_cnt: 0
}


begin

  #extract out to its own file?
  #Takes clade standards and makes a mixturized consensus standard.
  if(Settings['optimization-mixturize-subtypes'])
    s_genes = standards.map(){|e| e.id.split(' ')[2]}.uniq.sort()
    s_subtypes = standards.map(){|e| e.id.split(' ')[1]}.uniq.sort()
    new_standards = []
    s_genes.each do |gene|
      s_subtypes.each do |subtype|
        hxb2 = standards.find() {|e| sp = e.id.split(' ')
          sp[1] == subtype and sp[2] == gene and sp[3] == 'HXB2'
        }
        seqs = standards.find_all(){|e| sp = e.id.split(' ')
          sp[1] == subtype and sp[2] == gene and sp[3] == 'STD'
        }

        next if(!hxb2)

        mixturized_seq = ''
        0.upto(hxb2.nucleotides.size() - 1) do |i|
          nucs = seqs.map(){|e| NUC_MIX[e.nucleotides[i]]}.flatten().uniq().sort()
          mixturized_seq += NUC_MIX.invert[nucs]
        end

        new_standards << Sequence.new(id: "ref_#{subtype} #{subtype} #{gene} HXB2", nuc: hxb2.nucleotides, reference: true, source_filepath: hxb2.source_filepath)
        new_standards << Sequence.new(id: "ref_#{subtype} #{subtype} #{gene} STD", nuc: mixturized_seq, reference: false, source_filepath: hxb2.source_filepath)
      end
    end
    standards = new_standards
  end



  batch_files.each_with_index do |batch, batch_index|
    failure_state = false
    errors = []

    begin

      puts "Processing file #{batch_index + 1} of #{batch_files.size}:  #{batch[:input]}"
      sample = Sample.new()
      sample.id = batch[:input].gsub(/^.+\//, '')
      sample.source_filepath = batch[:input]

      #lets try playing around with a file.
      fastq = Fastq.new(batch[:input])

      #assign samplenames for batch output.
      if(options[:use_sampleids])
        invalid_characters = /[ \/:*?"<>|\\]/  #gets rid of possible bad filename characters for windows and unix.
        if(fastq.data.first())
          sampleid = fastq.data.first().annotations['sampleid'].gsub(invalid_characters, '')
          barcode = fastq.data.first().annotations['barcode'].gsub(invalid_characters, '')
          output = options[:batch_output] + '/' + sampleid + '+' + barcode
          if(batch_files.find(){|e| e[:output] == output })
            puts "Output name #{output} already exists, output label will use filename instead."
          else
            batch[:output] = output
          end
        else
          #use the filename previously assigned.
        end
      end

      job_add_num_sequences(job, fastq.data.size())
      puts "Loaded #{fastq.data.size()} sequences."

      #Gene setup
      gene_list.each do |gene_name|
        hxb2_standard = standards.find() do |e|
          e.id.split(' ')[2] == gene_name and e.id.split(' ')[3] == 'HXB2'
        end
        sample.genes << Gene.new(name: gene_name, hxb2_standard: hxb2_standard)
      end


      #prepping for alignments
      quit_gene_early = {}
      quit_gene_early.default = false
      gene_stds_hash = {}
      sample.genes.each do |sample_gene|
        gene_stds = standards.find_all() do |e|
          e.id.split(' ')[2] == sample_gene.name and e.id.split(' ')[3] == 'STD'
        end
        gene_stds_hash[sample_gene.name] = gene_stds
      end
      fastq_data = fastq.data
      if(options[:inspect])
        fastq_data = [fastq_data.find(){|fq| fq.id == options[:inspect]}]
      end
      if(Settings['debug'])
        fastq_data = fastq_data[0 .. Settings['debug-seq-limit']]
      end


      #Main alignment code
      fastq_data.each_with_index do |sequence, dex|
        sample.genes.each do |sample_gene|
          next if(quit_gene_early[sample_gene.name])
          gene_stds = gene_stds_hash[sample_gene.name]

          #start here
          alignments = [] #list of decent alignments generated
          seq_clipped = nil #for optimization if enabled

          gene_stds.each() do |standard|
            clade = standard.id.split(' ')[0]
            ref_standard = standard
            hxb2_standard = standards.find(){|e| e.id == standard.id.gsub(/STD$/, 'HXB2') }

            alignment = Alignment.new()
            if(Settings['optimization-clipped-alignment'] and seq_clipped)
              alignment.align(std: ref_standard, seq: seq_clipped)
            else
              alignment.align(std: ref_standard, seq: sequence)
            end
            match_perc = alignment.match_perc()

            if(match_perc < 0.10)
              alignment = Alignment.new()
              alignment.align(std: ref_standard, seq: sequence.reverse_complement())
              #puts "Bad alignment [#{sample_gene.name}] (#{match_perc}), r_complement: #{alignment.match_perc()}."
              match_perc = alignment.match_perc()
            end

            if(match_perc < Settings['optimization-stop-poor-alignments-threshold'] and
              Settings['optimization-stop-poor-alignments'])
              break #quick exit, as this doesn't align at all.
            #elsif(match_perc_parts < 0.45)
            #  puts "Poor partial match #{match_perc} #{match_perc_parts}"
            #  puts alignment.details()
            #  break
            elsif(match_perc > Settings['reference-match-threshold'])
              if(Settings['optimization-clipped-alignment'] and seq_clipped == nil)
                #I guess its fine because indels are preserved and the alignment would get rid of dashes anyway...
                tb = Settings['optimization-clipped-alignment-buffer']
                trim_start = ((alignment.trim_start - tb) > 0 ? (alignment.trim_start - tb) : 0)
                trim_end = ((alignment.trim_end + tb) < alignment.seq.nucleotides.size() ? (alignment.trim_end + tb) : alignment.seq.nucleotides.size() - 1)
                seq_clipped = alignment.seq.copy()
                seq_clipped.nucleotides = alignment.seq.nucleotides[trim_start .. trim_end]
              end

              #TODO, get rid of
              #puts "clade #{clade} match_perc #{match_perc} #{alignment.match_perc_parts(60)}"

              alignments << { ref_standard: ref_standard,
                clade: clade,
                alignment: alignment,
                match_perc: match_perc,
                hxb2_standard: hxb2_standard}
            else
              #puts "Bad I guess?  #{match_perc}"
            end

          end #end gene_stds.each

          #pick best alignment
          best_alignment = alignments.sort(){|a,b| b[:match_perc] <=> a[:match_perc]}.first()
          if(best_alignment == nil)
            stats[:no_match_cnt] += 1
            #puts "#{sequence.id}:  No Match"
          else
            stats[:reverse_cnt] += 1 if(best_alignment[:alignment].rev_comp)
            stats[:forward_cnt] += 1 if(!best_alignment[:alignment].rev_comp)
            stats[:total_cnt] += 1
            #puts "#{sequence.id}:  #{best_alignment[:clade]} #{(best_alignment[:match_perc] * 100).round(2)}%"
            #okay, put these alignments in the appropriate clade.
            subtype = best_alignment[:ref_standard].id.split(' ')[1]
            gene_clade = sample_gene.clades.find(){|e| e.id == best_alignment[:ref_standard].id}
            if(gene_clade.nil?)
              gene_clade = Clade.new()
              gene_clade.id = best_alignment[:ref_standard].id
              gene_clade.subtype = subtype
              gene_clade.ref_standard = best_alignment[:ref_standard]
              gene_clade.hxb2_standard = best_alignment[:hxb2_standard]
              sample_gene.clades << gene_clade
            end

            #re-align (#technically we only need to do it if our optimizations messed it up)
            final_alignment = nil
            if(Settings['optimization-clipped-alignment'])
              final_alignment = Alignment.new()
              if(best_alignment[:alignment].rev_comp)
                final_alignment.align(std: best_alignment[:ref_standard], seq: sequence.reverse_complement())
              else
                final_alignment.align(std: best_alignment[:ref_standard], seq: sequence)
              end
            else
              final_alignment = best_alignment.alignment
            end


            #alignment post processing here----------------------?

            if(Settings['reject-poor-aligned-ends']) #experimental, disabled by default
              good = reject_poor_aligned_ends(alignment: final_alignment)
              next if(!good) #skip poor alignments
            end

            if(Settings['reduce-dual-homopolymers']) #experimental, disabled by default
              reduce_dual_homopolymers(alignment: final_alignment )
            end

            if(Settings['optimization-coverage-limit'])
              if(optimization_coverage_limit(gene: sample_gene, alignment: final_alignment))
                quit_gene_early[sample_gene.name] = true
              end
            end

            #alignment ready
            gene_clade.alignments << final_alignment
          end

          #final sequence
          if(Settings['optimization-coverage-limit'] and
            fastq_data.size() - 1 == dex and
            !quit_gene_early[sample_gene.name])

            puts "#{sample_gene.name} - Target coverage not reached (min cov: #{sample_gene.aa_coverage.min} target: #{Settings['optimization-coverage-limit-target']})."
          end

        end #end sample.genes.each



        #Update job info every now and then.
        job_save(job) if(job_update(dex, job, sample))

        #quit if we are all done every gene.
        if(sample.genes.map(){|sg| quit_gene_early[sg.name] }.all?(){|a| a})
          break
        end

      end #end fastq.data.each


      #second pass
      sample.genes.each do |sample_gene|
        #trim gene edges sequence improvement.
        trim_gene_edges(gene: sample_gene) if(Settings['trim-gene-edges'])

        #homopolymer fixes/realignment improvements
        homopolymer_fixes(gene: sample_gene) if(Settings['homopolymer-fixes'])

        if(Settings['alignment-optimize'])
          ts = Time.now
          alignment_optimize(gene: sample_gene)
          #puts "align-optimize took #{(Time.now - ts)} seconds"
        end


        if(Settings['sweep-inserts']) #experimental, disabled by default
          sample_gene.clades.each do |clade|
            clade.alignments.each do |final_alignment|
              sweep_inserts(alignment: final_alignment)
            end
          end
        end

      end


      #clean this up or move elsewhere
      if(options[:inspect])
        File.open(batch[:output] + '.inspect.txt', 'w') do |file|
          sample.genes.each do |gene|
            gene.subtypes.each do |clade|
              #write out trimmed raw alignments
              alignment = clade.alignments.first()

              qual_array = []
              alignment.trim_start.upto(alignment.trim_end) do |i|
                qual_array << alignment.qual_at(pos: i)
              end
              qual_str = qual_to_string(qual_array)

              file.puts ">#{gene.name} #{clade.id} #{alignment.std.id}"
              file.puts alignment.std.nucleotides[alignment.trim_start  .. alignment.trim_end]
              file.puts ">#{alignment.seq.id}"
              file.puts alignment.seq.nucleotides[alignment.trim_start  .. alignment.trim_end]
              file.puts qualstr
              file.puts

              str = ''
              str += "Insertions: " + alignment.insertions.map(){|e|
                if(e.size() == 1)
                  alignment.loc_ref(e.first).to_s + ':' + alignment.seq.nucleotides[e]
                else
                  alignment.loc_ref(e.first).to_s + '-' + (alignment.loc_ref(e.first) + e.size()).to_s + ':' + alignment.seq.nucleotides[e]
                end
              }.join(' ') + "\n"

              str += "Deletions: " + alignment.deletions.map(){|e|
                if(e.size() == 1)
                  alignment.loc_ref(e.first).to_s
                else
                  alignment.loc_ref(e.first).to_s + '-' + (alignment.loc_ref(e.first) + e.size()).to_s
                end
              }.join(' ') + "\n\n"
              file.puts str
            end
          end
        end

      end


      #We should choose to accept or fail it here I think
      #If we fail, we don't need reports.

      #1:  Did we reach the target count for each gene?
      #What if we didn't bother with integrase or something?  Hrm.

      sample.genes.each do |sample_gene|
        if(sample_gene.aa_coverage.min < Settings['optimization-coverage-limit-target'])
          errors << "Did not reach the target coverage for the #{sample_gene.name} region."
          failure_state = true
        end
      end

  #    pp stats
      if(Settings['debug'])
        sample.print_stats()
        #puts report.text_report
      end

      if(Settings['debug']) #fasta of less-nice insertion sequences.
        fastasave = []
        sample.genes.each do |gene|
          gene.clades.each do |clade|
            clade.alignments.each do |alignment|
              if(alignment.insertions.size() > 0)
                fastasave << [alignment.std.id, alignment.std.nucleotides[alignment.trim_start  .. alignment.trim_end]]
                fastasave << [alignment.seq.id, alignment.seq.nucleotides[alignment.trim_start  .. alignment.trim_end]]
              end
            end
          end
        end
        save_fasta(filename: batch[:output] + ".debug_ins.fas", sequences: fastasave)
      end

    rescue
      failure_state = true
      errors << "Processing error:  #{$!.to_s}"
      puts $!
      puts $!.backtrace
    end



    #Save job here instead?
    if(failure_state == false)
      #Output reports here

      report = ResistanceReport.new(label: sample.id)
      report_for = ResistanceReport.new(label: sample.id)
      report_rev = ResistanceReport.new(label: sample.id)

      sample.genes.each do |gene|
        filename_root = batch[:output] + "." + gene.name
        gene_alignments = []
        gene.clades.each {|clade| gene_alignments += clade.alignments}
        gene_alignments_for = gene_alignments.select(){|a| a.rev_comp == false}
        gene_alignments_rev = gene_alignments.select(){|a| a.rev_comp == true}

        if(Settings['debug'])
          MutationFreq.save_nuc(
            filename: filename_root + ".nuc-freq.csv",
            alignments: gene_alignments,
            hxb2_standard: gene.hxb2_standard,
            trim_range: gene.trim_range)
        end
        MutationFreq.save_aa(
          filename: filename_root + ".aa-freq.csv",
          alignments: gene_alignments,
          hxb2_standard: gene.hxb2_standard,
          trim_range: gene.trim_range )
        if(Settings['output-filter-direction'])
          #filter by direction.
          MutationFreq.save_aa(
            filename: filename_root + ".aa-freq.for.csv",
            alignments: gene_alignments_for,
            hxb2_standard: gene.hxb2_standard,
            trim_range: gene.trim_range )
          MutationFreq.save_aa(
            filename: filename_root + ".aa-freq.rev.csv",
            alignments: gene_alignments_rev,
            hxb2_standard: gene.hxb2_standard,
            trim_range: gene.trim_range )
        end

        #insertion lists
        File.open(filename_root + ".insertions.csv", 'w') do |file|
          data = MutationFreq.get_ins_frequencies(alignments: gene_alignments)
          file.puts "NUC_LOC,AA_LOC,CNT,PERC,SIZE,NUC,AA,VALID"
          data.each do |dat|
            if(dat[2].to_f / gene_alignments.size().to_f > 0.10) #10% threshold
              file.puts dat.join(',')
            end
          end
        end

        #clean fastas
        fastasave = []
        fastasave_for = []
        fastasave_rev = []
        fastasave << ['HXB2', gene.hxb2_standard.nucleotides]
        fastasave_for << ['HXB2', gene.hxb2_standard.nucleotides]
        fastasave_rev << ['HXB2', gene.hxb2_standard.nucleotides]

        gene_alignments.each do |alignment|
          fastasave << [alignment.seq.id, alignment.seq_clean[gene.trim_range]]
          if(alignment.rev_comp)
            fastasave_rev << [alignment.seq.id, alignment.seq_clean[gene.trim_range]]
          else
            fastasave_for << [alignment.seq.id, alignment.seq_clean[gene.trim_range]]
          end
        end
        save_fasta(filename: filename_root + ".fas", sequences: fastasave)
        if(Settings['output-filter-direction'])
          save_fasta(filename: filename_root + ".for.fas", sequences: fastasave_for)
          save_fasta(filename: filename_root + ".rev.fas", sequences: fastasave_rev)
        end

        #Add data to resistance report
        aa_consensus = MutationFreq.get_consensus_aa_array(
          alignments: gene_alignments, trim_range: gene.trim_range)
        aa_consensus_for = MutationFreq.get_consensus_aa_array(
          alignments: gene_alignments_for, trim_range: gene.trim_range)
        aa_consensus_rev = MutationFreq.get_consensus_aa_array(
          alignments: gene_alignments_rev, trim_range: gene.trim_range)

        report.add_aa(aa_seq: aa_consensus, gene: gene.name,
          aa_hxb2: nuc_to_aa_array(gene.hxb2_standard.nucleotides) )
        report_for.add_aa(aa_seq: aa_consensus_for, gene: gene.name,
          aa_hxb2: nuc_to_aa_array(gene.hxb2_standard.nucleotides) )
        report_rev.add_aa(aa_seq: aa_consensus_rev, gene: gene.name,
          aa_hxb2: nuc_to_aa_array(gene.hxb2_standard.nucleotides) )

        #nucleotide consensus
        nuc_consensus = MutationFreq.get_consensus_nuc(
          alignments: gene_alignments, trim_range: gene.trim_range )
        File.open(filename_root + ".consensus.txt",'w') do |file|
          file.puts nuc_consensus
        end

      end

      #Save resistance report
      File.open(batch[:output] + ".resistance-report.txt", 'w') do |file|
        file.puts report.text_report
      end
      if(Settings['output-filter-direction'])
        File.open(batch[:output] + ".resistance-report.for.txt", 'w') do |file|
          file.puts report_for.text_report
        end
        File.open(batch[:output] + ".resistance-report.rev.txt", 'w') do |file|
          file.puts report_rev.text_report
        end
      end

      #save job
      job_update(fastq.data.size() - 1, job, sample, true) #force final update
      job_complete(job, report)  #add final report
      job_save(job)
      job_zip(job, batch[:output]) if(options[:zip])

      puts "Processing successful."
    else #failure_state==true
      #save error information
      job_error(job, errors.join("\n"))
      job_save(job)
      job_zip(job, batch[:output]) if(options[:zip])

      puts "Processing failed:"
      puts errors
    end


  end #end processing files/folders

rescue
  puts $!
  puts $!.backtrace
end


exit()
