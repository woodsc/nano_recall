=begin
lib/alg/aligner.rb
Copyright (c) 2007 BCCFE

Aligns a sequence (with an algorithm no less)

Maybe some of these methods should be elsewhere?  (Smelt and trim?)
=end

require 'lib/recall_config'
require 'lib/alg/_alignment' #c-extension
require 'lib/conversions'

class Aligner
  extend SeqConversions
	def Aligner.align(data, isweb=false)
		#Align each sequence to a standard(make a standard/sequence pair)
		
		aligned_data = [] # [[standard, seq], ...]
		
    reject_alignment_percent = RecallConfig['aligner.reject_alignment_percent'].to_f
    blankout_bad_edges = (RecallConfig['aligner.blankout_bad_edges'] == 'true')
        
		data.primers.each do |p|
      if(p.direction == nil)
        data.errors.push("Primer has unknown direction: #{isweb ? p.name : p.primerid}")
      end
            
      if(data.long_standard != nil)
        aligned_data.push([data.long_standard,p.called])
      else
        aligned_data.push([data.standard,p.called])
      end
		end

		#align each element of aligned_data
		aligned_data.each do |elem|
			Aligner.run_alignment(elem)
		end
        
    #Should probably trim hereish.
    trim_standard(aligned_data, data.standard) if(aligned_data.size != 0 and data.long_standard != nil)
        
       
#        puts aligned_data[0][0].join()
        
    
=begin
#This works better than the version below, but Richard doesn't really want changes to basecalling
#Must fail things that don't really match up
    aligned_data.each_with_index do |elem, dex|
      score = 0
      cnt = 0
      fb = nil
      lb = nil

      0.upto(elem[0].size - 1) do |i|
        fb = i if(fb == nil and elem[0][i] != '-' and elem[1][i] != '-')
        lb = i if(elem[0][i] != '-' and elem[1][i] != '-')
      end
      
      fb = 0 if(fb == nil)
      lb = (elem[0].size - 1) if(lb == nil)
      
      fb.upto(lb) do |i|
        if(elem[0][i] != '-' and elem[1][i] != '-' and elem[0][i] == elem[1][i])
          score += 1
          cnt += 1
        else#if(elem[1][i] != '-')
          cnt += 1
        end
      end
#      puts elem[0].join('')
#      puts elem[1].join('')
#      puts "#{data.primers[dex].name()}:\t#{score.to_f / cnt.to_f}"
#      File.open('primer_alignment.log', 'a') do |file|
#        file.puts "#{data.project} - #{data.sample} - #{data.primers[dex].primerid()}:\t#{score.to_f / cnt.to_f}"
#      end
      
#            puts (score.to_f / cnt.to_f)
      if(score.to_f / cnt.to_f < reject_alignment_percent or cnt < 15)
        elem[0] = nil #Meaning this was probably a misalignment
        elem[1] = nil
        data.errors.push("Failing primer #{isweb ? data.primers[dex].name : data.primers[dex].primerid}; can't align to standard")
        data.primers[dex] = nil
      end
    end
=end

    
    #Must fail things that don't really match up
    aligned_data.each_with_index do |elem, dex|
        score = 0
        cnt = 0
        0.upto(elem[0].size - 1) do |i|
          if(elem[0][i] != '-' and elem[1][i] != '-' and elem[1][i] != 'N'   and (@@ambig_nucs[elem[0][i]].all?{|b| @@ambig_nucs[elem[1][i]].include?(b)} or @@ambig_nucs[elem[1][i]].all?{|b| @@ambig_nucs[elem[0][i]].include?(b)}))
          #if(elem[0][i] != '-' and elem[1][i] != '-'  and @@ambig_nucs[elem[0][i,1]].all?{|b| @@ambig_nucs[elem[1][i,1]].include?(b)} or @@ambig_nucs[elem[1][i,1]].all?{|b| @@ambig_nucs[elem[0][i,1]].include?(b)})
          #if(elem[0][i] != '-' and elem[1][i] != '-' and elem[1][i] != 'N'   and elem[0][i] == elem[1][i])
                score += 1
                cnt += 1
            elsif(elem[0][i] != '-' and elem[1][i] != '-')
                cnt += 1
            end
        end
        
#            puts (score.to_f / cnt.to_f)
        if(score.to_f / cnt.to_f < reject_alignment_percent or cnt < 15)
            elem[0] = nil #Meaning this was probably a misalignment
            elem[1] = nil
            data.errors.push("Failing primer #{isweb ? data.primers[dex].name : data.primers[dex].primerid}; can't align to standard")
            data.primers[dex] = nil
        end
    end



    data.primers.delete_if {|p| p == nil }
    aligned_data.delete_if {|elem| elem[0] == nil}
    
        
		final = []
		#Align standards to each other(along with the sequence). Cool!
        
		run_alignment_merge(aligned_data) if(aligned_data.size != 0)
#		fine_align(aligned_data) if(aligned_data.size != 0)
    correct_alignment(aligned_data) if(aligned_data.size != 0)

    ProjectMethods.adjust_alignment(aligned_data) #Project custom code
    
    #get rid of stop codon matcher thingy
    #aligned_data.each do |elem|
    #    (elem[0].size - 1).downto(0) do |i|
    #        if(elem[0][i] == '$')
    #            elem[0].delete_at(i)
    #            elem[1].delete_at(i)
    #        end
    #    end
    #end
    
    

    data.standard = aligned_data[0][0] if(aligned_data.size != 0)

    i=0
		data.primers.each do |p|
      if(aligned_data[i][0] == nil)
          p.edit = nil
          i += 1
          next
      end
      
			p.edit = aligned_data[i][1] 

			#align the stuff to the edit(annoying)
			p.edit.each_with_index do |v, j|
				if(v == '-')
					p.called.insert(j, '-')
					p.uncalled.insert(j, '-')
					p.called_area.insert(j, '-')
					p.uncalled_area.insert(j, '-')
					p.loc.insert(j, '-')
					p.qual.insert(j, '-')
					p.ignore.insert(j, '-')
          p.amp_a.insert(j, '-')
          p.amp_c.insert(j, '-')
          p.amp_g.insert(j, '-')
          p.amp_t.insert(j, '-')
				end
			end
			
			i += 1
		end

    #Trim primer ends if it looks misaligned.
    if(blankout_bad_edges)
      data.primers.each do |p|
        #find primer start
        st = p.primer_start
        #If there are dashes in the first 10 bases, blank it out
        lastdash = p.edit[st, 10].rindex('-')
        if(lastdash)
          st.upto(st + lastdash) do |i|
            if(p.edit[i] != '-')
              #puts p.edit[i]
              p.ignore[i] = 'L'
            end
          end
        end
        #find primer end
        en = p.primer_end
        firstdash = p.edit[en - 10, 10].index('-')
        if(firstdash)
          ((en - 10) + firstdash).upto(en) do |i|
            if(p.edit[i] != '-')
              #puts p.edit[i]
              p.ignore[i] = 'L'
            end
          end
        end
      end
    end
	end
	
	
	#Runs the c++ alignment algorithm(might replace in ruby later for fun)
	def Aligner.run_alignment(elem)
		gap_init_penalty = RecallConfig['aligner.gap_initialization_penalty'].to_i
		gap_extend_penalty = RecallConfig['aligner.gap_extension_penalty'].to_i
		align_prog = RecallConfig['aligner.alignment_program']
		tmp = align_it(elem[0].join(''), elem[1].join(''), gap_init_penalty, gap_extend_penalty)
    elem[0] = tmp[0].split('')
    elem[1] = tmp[1].split('')
	end
	
  def Aligner.trim_stop_codon(data)
    #get rid of stop codon matcher thingy
    (data.standard.size - 1).downto(0) do |i|
      if(data.standard[i] == '$')
        data.standard[i] = '-'
        data.assembled[i] = '-'
=begin
      data.primers.each do |p|
          p.called[i] = '-'
          p.uncalled[i] = '-'
          p.called_area[i] = '-'
          p.uncalled_area[i] = '-'
          p.loc[i] = '-'
          p.qual[i] = '-'
          p.ignore[i] = '-'
          p.amp_a[i] = '-'
          p.amp_c[i] = '-'
          p.amp_g[i] = '-'
          p.amp_t[i] = '-'
      end
=end
      end
    end
  end
    
  #Need to write this code?
  def Aligner.trim_standard(aligned_data, standard)
      aligned_data.each do |elem|
          
          t_aligned = align_it(elem[0].join(''), standard.join(''),3,1)
          #t_aligned[1] =~ /^(-+)[^-]+(-+)$/
          t_aligned[1] =~ /^(-*)[^-]+(-*)$/
          
          #puts t_aligned[0]
          #puts t_aligned[1]
          #puts

          count_start = $1.size
          count_end = $2.size
      #kill everything at the start
          i = 0
          while(count_start != 0)
              if(elem[0][i] != '-')
                  count_start -= 1 
                  elem[0][i] = '-'
              end
              i+=1
          end
          i = 0
      #kill everything at the end
          while(count_end != 0)
              if(elem[0][-i] != '-')
                  count_end -= 1 
                  elem[0][-i] = '-'
              end
              i+=1
          end
      end
  end
    
	#Note, nothing before or after the standard is 
	#aligned.  (it may turn out that way, but its a coincidence)
	def Aligner.run_alignment_merge(aligned_data)
		1.upto(aligned_data.size - 1) do |i|
			#align 0 to i-1 to i.  0 to i - 1 should already be aligned.
			j = 0
 
			#Quicky alignment to align the start of standards
			#This is mostly so the unaligned parts look a little more aligned
			aligned_data[0][0].join('') =~ /^(\-*)/
			sizea = $1.length
			aligned_data[i][0].join('') =~ /^(\-*)/
			sizeb = $1.length

			if(sizea > sizeb)
				aligned_data[i][0].insert(0, *(['-'] * (sizea - sizeb)))
				aligned_data[i][1].insert(0, *(['-'] * (sizea - sizeb)))
			elsif(sizea < sizeb)
				aligned_data[0 .. i -1].each do |elem|
					elem[0].insert(0, *(['-'] * (sizeb - sizea)))
					elem[1].insert(0, *(['-'] * (sizeb - sizea)))
				end
			end
			
			#Ok, on with the main alignment
			while(j < aligned_data[0][0].size or j < aligned_data[i][0].size)
				if(aligned_data[0][0][j] == aligned_data[i][0][j])
					#OK
					j += 1
				elsif(aligned_data[0][0][j] == '-' and aligned_data[i][0][j] != '-')
					#Ok, add a dash to the aligned_data[i]
					aligned_data[i][0].insert(j, '-')
					aligned_data[i][1].insert(j, '-')
				elsif(aligned_data[0][0][j] != '-' and aligned_data[i][0][j] == '-')
					#Ok, add a dash to the aligned_data[0 .. i - 1]
					aligned_data[0 .. i -1].each do |elem|
						elem[0].insert(j, '-')
						elem[1].insert(j, '-')
					end
				else
					#This shouldn't happen.
					puts "Standards don't match?  Gameover man! GAMEOVER!"
				end
				
			end
		end
	end
=begin	
	#No longer used
	def Aligner.fine_align(data)
		2.upto(data[0][0].length - 3) do |i|
			data.each do |p|
				if(p[1][i] == '-' and p[1][i - 1] != '-' and p[1][i + 1] != '-')
					if(p[0][i - 1] == '-' and p[0][i - 2] != '-' and p[0][i] != '-')
						#Shift dash backwards
						p[1][i] = p[1][i - 1]
						p[1][i - 1] = '-'
						puts "Fine aligned"
					elsif(p[0][i + 1] == '-' and p[0][i] != '-' and p[0][i + 2] != '-')
						#Shift dash forwards
						p[1][i] = p[1][i + 1]
						p[1][i + 1] = '-'
						puts "Fine aligned"
					elsif(p[0][i - 2] == '-' and p[0][i - 3] != '-' and p[0][i - 1] != '-' and p[1][i - 2] != '-')
						#Shift dash backwards 2
						p[1][i] = p[1][i - 1]
						p[1][i - 1] = p[1][i - 2]
						p[1][i - 2] = '-'
						puts "Fine aligned x 2"
					elsif(p[0][i + 2] == '-' and p[0][i + 1] != '-' and p[0][i + 3] != '-' and p[1][i + 2] != '-')
						#Shift dash forwards 2
						p[1][i] = p[1][i + 1]
						p[1][i + 1] = p[1][i + 2]
						p[1][i + 2] = '-'
						puts "Fine aligned x 2"
					end
				end
			end
		end
	end
=end    
    #This might be a better way to do the fine_align
    def Aligner.correct_alignment(data)
        #This should find single base insertions and line them up to single base deletions.
        
        
        1.upto(data[0][0].length - 2) do |i|
            if(data[0][0][i] == '-' and data[0][0][i - 1] != '-' and data[0][0][i + 1] != '-')
                #initiate search around the dash.
                #Look around for 15 bases.
                pos = nil
                n = nil
            
                (i - 1).downto(i - 15) do |j|
                    if(data[0][0][j] != '-' and j > -1 and pos == nil)
                        n = data.find_all do |d|
                            if(d[1][j] == '-' and d[1][j - 1] != '-'  and d[1][j + 1] != '-')
                                true
                            else
                                false
                            end
                        end
                        pos = j if(n.size >= 2)
                    end
                end
                
                (i + 1).upto(i + 15) do |j|
                    if(data[0][0][j] != '-' and j < data[0][0].length and pos == nil)
                        n = data.find_all do |d|
                            if(d[1][j] == '-' and d[1][j - 1] != '-'  and d[1][j + 1] != '-')
                                true
                            else
                                false
                            end
                        end
                        pos = j if(n.size >= 2)
                    end
                end

                #Now we want to shift them all
                if(pos)
#                    puts data[0][0][i - 15 .. i + 15].join('')
#                    puts '~' * 31
                    n.each do |d|
#                        puts d[1][i - 15 .. i + 15].join('')
                        #Must shift pos to i
                        if(pos > i)
                            d[1][i + 1 .. pos] = d[1][i .. pos - 1]
                            d[1][i] = '-'
                        elsif(pos < i)
                            d[1][pos .. i - 1] = d[1][pos + 1 .. i]
                            d[1][i] = '-'
                        end
#                        puts d[1][i - 15 .. i + 15].join('')
#                        puts
                    end
                end
                
            end
        end
        
        #An insert is when all(most?) of the edits have a letter, but standard has a dash
        
    end

end

