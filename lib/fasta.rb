require "#{File.dirname(__FILE__)}/sequence.rb"
require "#{File.dirname(__FILE__)}/settings.rb"

def load_fasta(filename: , reference: false)
  sequences = []

  File.open(filename, 'r') do |file|
    label = nil
    seq = nil
    file.each_line do |line|
      if(line =~ /^>/)
        if(label != nil)
          sequences << Sequence.new(id: label, nuc: seq, reference: reference, source_filepath: filename)
        end
        label = line[1 .. -1].strip()
        seq = ''
      else
        seq += line.strip()
      end
    end
    if(label != nil and seq != nil and seq != '')
      sequences << Sequence.new(id: label, nuc: seq, reference: reference, source_filepath: filename)
    end
  end

  return sequences
end

def save_fasta(filename: , sequences: )
  File.open(filename, 'w') do |file|
    sequences.each do |seq|
      if(seq.class == Sequence)
        file.puts ">#{seq.label}"
        file.puts "#{seq.nucleotides}"
      elsif(seq.class == Array)
        file.puts ">#{seq[0]}"
        file.puts "#{seq[1]}"
      end
    end
  end
end
