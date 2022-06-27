require "../../lib/alignment/_alignment.rb"

final = [] #[ref_name, region, hxb2_align, sequence_align]
hxb2_ref = nil
raw_standards = []

hxb2_genes = { #[start, end]
  :protease => [2253,2549],
  :rt => [2550, 3869],
  :int => [4230, 5093],
}

puts "prepping"

#load hiv.lanl.gov reference standards.
File.open("HIV1_REF_2019_genome_DNA.csv") do |file|
  file.each_line do |line|
    raw_standards << line.strip.split(',')
  end
end

#load HXB2 sequence.
File.open("HIV-HXB2.txt") do |file|
  hxb2_ref = file.gets(nil).strip()
end

#filter to subtypes we care about
target_standards = raw_standards.filter() do |e|
  if(e[0] =~ /Ref[\.\_](A\d)[\.\_]/)
    true
  elsif(e[0] =~ /Ref[\.\_](B)[\.\_]/)
    true
  elsif(e[0] =~ /Ref[\.\_](C)[\.\_]/)
    true
  elsif(e[0] =~ /Ref[\.\_](D)[\.\_]/)
    true
  elsif(e[0] =~ /Ref\.01[\.\_](AE)[\.\_]/)
    true
  elsif(e[0] =~ /Ref\.02[\.\_](AG)[\.\_]/)
    true
  else
    false
  end
end

#final = [] #[ref_name, region, hxb2_align, sequence_align]
#hxb2_genes
puts "aligning" #slowest bit
target_standards.each do |ref|
  result = align_it(hxb2_ref, ref[1], 3, 1)
  ref << result #ref[2] contains the alignment
end

puts "cutting"
#cut them up
target_standards.each do |ref|
  hxb2_dex = 0
  gene_cuts = {}

  0.upto(ref[2][0].size() - 1) do |i|
    hxb2_genes.each do |gene, gene_range|
      gene_cuts[gene] = ['',''] if(i == 0) #init gene_cuts

      if((gene_range[0] - 1) <= hxb2_dex and (gene_range[1] - 1) >= hxb2_dex)
        gene_cuts[gene][0] += ref[2][0][i]
        gene_cuts[gene][1] += ref[2][1][i]
      end

    end

    hxb2_dex += 1 if(ref[2][0][i] != '-')
  end

  gene_cuts.each do |gene, gene_cut|
    subtype = 'X'
    if(ref[0] =~ /Ref[\.\_](A\d)[\.\_]/)
      subtype = $1
    elsif(ref[0] =~ /Ref[\.\_](B)[\.\_]/)
      subtype = $1
    elsif(ref[0] =~ /Ref[\.\_](C)[\.\_]/)
      subtype = $1
    elsif(ref[0] =~ /Ref[\.\_](D)[\.\_]/)
      subtype = $1
    elsif(ref[0] =~ /Ref\.01[\.\_](AE)[\.\_]/)
      subtype = $1
    elsif(ref[0] =~ /Ref\.02[\.\_](AG)[\.\_]/)
      subtype = $1
    end
    final << [ref[0], gene, gene_cut[0], gene_cut[1], subtype]
  end

end

puts "saving"

File.open("reference_standards.fas", 'w') do |file|
  final.each do |elem|
    file.puts ">#{elem[0]} #{elem[4]} #{elem[1]} HXB2\n#{elem[2]}"
    file.puts ">#{elem[0]} #{elem[4]} #{elem[1]} STD\n#{elem[3]}"
  end
end


#lets see if there is any point in keeping all these triplicates.
#Result:  Turns out there IS more differences than I would have expected, so lets keep them for now.
=begin
hxb2_genes.each do |gene, range|
  puts "Gene: #{gene}"

  final.find_all(){|e| e[1] == gene }.each do |fin|
    clade = fin[0].split('.')[1]
    #compare similar clades
    final.find_all(){|e| e[1] == gene and e[0] != fin[0] and e[0].split('.')[1] == clade }.each do |fin2|
      diff = 0
      0.upto(fin[3].size() - 1) do |i|
        diff += 1 if(fin[3][i] != fin2[3][i])
      end

      puts "#{fin[0]} - #{fin2[0]}:  #{diff}"
    end
  end
end
=end

#
