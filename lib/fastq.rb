#fastq.rb converts fastq files into sequence objects.
#
#sample:
=begin
@b6f3c08d-4848-44db-8d3d-98ecbfea3378 runid=af5ea6be5adb2f9f887925079cffeb55c523
804b read=33 ch=56 start_time=2020-02-06T07:32:42Z flow_cell_id=ABY287 protocol_
group_id=HIV sample_id=MIN00021 barcode=barcode01
TGAGTACGCTTCGTTCAGTTA....BLAH
+
#&#%$$%%,5786D5:3700-+*0//31:8893B:7'')(,12%/6:8,,,11.,*///00,,(%((--35+'2
=end

require 'zlib'
require 'stringio'
require "#{File.dirname(__FILE__)}/sequence.rb"
require "#{File.dirname(__FILE__)}/settings.rb"

class Fastq
  attr_accessor :data

  def initialize(filename)
    @data = []
    sequence = Sequence.new()
    sequence.source_filepath = filename

    file_lines = []

    #A bit complicated because rubys gzip library is dumb as balls.
    File.open(filename, 'rb') do |file|
      if(filename =~ /\.gz$/)
        file_str = file.read()
        io = StringIO.new(file_str)
        done = false
        while(!done)
          gz = Zlib::GzipReader.new(io)
          gz.each_line() do |line|
            file_lines.push(line)
          end
          if(gz.unused)
            io.pos -= gz.unused.size()
          elsif(io.eof)
            done = true
          end
          gz.finish()
        end
      else
        file.each do |line|
          file_lines << line
        end
      end
    end

    #puts "Lines: " + file_lines.size().to_s
    #exit()

    file_lines.each_with_index do |line, line_num|
      if(line_num % 4 == 0 && line.strip() == '')
        #pass, empty line, probably at the end.
      elsif(line_num % 4 == 0)
        raise "Invalid fastq format line #{line_num}" if(line[0,1] != '@')
        tmp = line.strip.split(' ')
        sequence.id = tmp[0][1 .. -1]
        tmp[1 .. -1].each() do |e|
          sequence.annotations[e.split('=')[0]] = e.split('=')[1]
        end
      elsif(line_num % 4 == 1)
        sequence.nucleotides = line.strip()
      elsif(line_num % 4 == 2)
        raise "Invalid fastq format line #{line_num}" if(line[0,1] != '+')
      elsif(line_num % 4 == 3)
        sequence.quality = line.strip().bytes.map(){|a| a - 33}
        @data.push(sequence)

        sequence = Sequence.new()
        sequence.source_filepath = filename
      end

      line_num += 1
    end

  end

end

#test
=begin
require 'pp'
fastq = Fastq.new("/home/woodsc/Downloads/ABY287_pass_barcode01_af5ea6be_0.fastq")
pp fastq.data[8].annotations
puts fastq.data.size()
puts fastq.data[8].nucleotides.size()
puts fastq.data[8].quality.size()
=end
