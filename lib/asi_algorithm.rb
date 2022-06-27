=begin
This class loads an ASI2 algorithm XML file and uses the loaded algorithm to
interpret amino acid strings into resistance scores.
=end


require "rexml/document"
#require 'jcode'


#code proper
class BNFVal
  attr_accessor :truth, :cond, :score, :logic, :flags
  def initialize(cond, truth=false, score=0, flags=[])
    @truth = truth
    @cond = cond
    @score = score
    @flags = flags
    @logic = nil #will be 'AND' or 'OR' if its from a condition2
  end
end

class AsiAlgorithm
  attr_accessor :debug, :alg_version, :alg_name, :level_def

  def initialize() #disabled algorithm, so I can use the method calls.

  end

  #slow, because REXML sucks.
  def initialize(filename, basic = false) #filename is the xml file containing the asi algorithm.  If null, use the latest HIVDB version.
    if(filename == nil)
      list = Dir[File.expand_path(File.dirname(__FILE__)) + "/HIVDB_*.xml"]
      list.sort!() #Sort and get the most recent.
      filename = list.last()
	  end

    if(!basic)
  		@pr_std = 'PQVTLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF'
  		@rt_std = 'PISPIETVPVKLKPGMDGPKVKQWPLTEEKIKALVEICTEMEKEGKISKIGPENPYNTPVFAIKKKDSTKWRKLVDFRELNKRTQDFWEVQLGIPHPAGLKKKKSVTVLDVGDAYFSVPLDEDFRKYTAFTIPSINNETPGIRYQYNVLPQGWKGSPAIFQSSMTKILEPFRKQNPDIVIYQYMDDLYVGSDLEIGQHRTKIEELRQHLLRWGLTTPDKKHQKEPPFLWMGYELHPDKWTVQPIVLPEKDSWTVNDIQKLVGKLNWASQIYPGIKVRQLCKLLRGTKALTEVIPLTEEAELELAENREILKEPVHGVYYDPSKDLIAEIQKQGQGQWTYQIYQEPFKNLKTGKYARMRGAHTNDVKQLTEAVQKITTESIVIWGKTPKFKLPIQKETWET'
  		@int_std = 'FLDGIDKAQDEHEKYHSNWRAMASDFNLPPVVAKEIVASCDKCQLKGEAMHGQVDCSPGIWQLDCTHLEGKVILVAVHVASGYIEAEVIPAETGQETAYFLLKLAGRWPVKTIHTDNGSNFTGATVRAACWWAGIKQEFGIPYNPQSQGVVESMNKELKKIIGQVRDQAEHLKTAVQMAVFIHNFKRKGGIGGYSAGERIVDIIATDIQTKELQKQITKIQNFRVYYRDSRNPLWKGPAKLLWKGEGAVVIQDNSDIKVVPRRKAKIIRDYGKQMAGDDCVASRQDED'

  		#Algorithm info
  		@debug = false
  		@alg_version = ''
  		@alg_name = ''

  		#definitions
  		@gene_def = [] #[['PR', ['PI']], ['RT',['NNRTI','NRTI']], ...]
  		@level_def = [] #[['1', 'Susceptible', 'S'], ...]
  		@drug_class = [] #[ ['PI', ['FPV/r', 'IDV/r', ...] ], ...]
  		@global_range = []  #[ ['-INF', '9', '1'] , ...]  #first two are the range, the third one is the res level
  		@comment_def = [] #something...

  		@drugs = []  #sub objects?  Hashes?
  		@mutation_comments = []  #maybe skip for now?  We don't really use this atm.
  		@mutations = []   #only for hivdb.  This is technically a hack that isn't part of the alg

  		#parse the file
  		dom = nil
  		File.open(filename) do |file|
  		  dom = REXML::Document.new(file)
  		end

  		#algorithm info
  		@alg_name = dom.elements["*/ALGNAME"].text()
  		@alg_version = dom.elements["*/ALGVERSION"].text()

  		#definitions
  		defs = dom.elements["*/DEFINITIONS"]

  		defs.elements.each("GENE_DEFINITION") do |elem|
  		  a = elem.elements['NAME'].text().strip()
  		  b = elem.elements['DRUGCLASSLIST'].text().split(',').map{|e| e.strip()}
  		  @gene_def << [a, b]
  		end
  		#puts @gene_def.inspect

  		defs.elements.each("LEVEL_DEFINITION") do |elem|
  		  a = elem.elements['ORDER'].text().strip()
  		  b = elem.elements['ORIGINAL'].text().strip()
  		  c = elem.elements['SIR'].text().strip()
  		  @level_def << [a, b, c]
  		end
  		#puts @level_def.inspect

  		defs.elements.each("DRUGCLASS") do |elem|
  		  a = elem.elements['NAME'].text().strip()
  		  b = elem.elements['DRUGLIST'].text().split(',').map{|e| e.strip()}
  		  @drug_class << [a, b]
  		end
  		#puts @drug_class.inspect

  		if(defs.elements["GLOBALRANGE"] != nil)
  		  @global_range = defs.elements['GLOBALRANGE'].text().strip().split(',').map{|e| e.gsub(/^[ \(\)\n]*/,'').gsub(/[ \(\)\n]*$/, '') }
  		  @global_range.map! do |e|
  			e =~ /\s*(\S+)\s*TO\s*(\S+)\s*=>\s*(\S+)\s*/
  			[$1,$2,$3]
  		  end
  		end
  		#puts @global_range.inspect

  		if(defs.elements["COMMENT_DEFINITIONS"] != nil)
  		  defs.elements.each("COMMENT_DEFINITIONS/COMMENT_STRING") do |elem|
  			a = elem.attributes.find{|a| a[0].upcase == 'ID'}[1].strip()
  			b = elem.elements['TEXT'].text().strip()
  			c = elem.elements['SORT_TAG'].text().strip()
  			@comment_def << [a,b,c]
  		  end
  		end
  		#puts @comment_def.inspect

  		#DRUGS
  		dom.elements.each("*/DRUG") do |drug_node|
  		  name = drug_node.elements['NAME'].text().strip()
  		  fullname = ''
  		  drug_node.elements.each("FULLNAME") do |elem|
  			fullname = elem.text().strip() #sometimes doesn't exist
  		  end

  		  rules = []
  		  drug_node.elements.each("RULE") do |rule_node|
  			condition = rule_node.elements['CONDITION'].text().strip().gsub(/\s+/, ' ')
  			actions = []

  			rule_node.elements.each("ACTIONS") do |action_node|
  			  action_node.elements.each("LEVEL") do |level|
  				actions << ['level', level.text().to_i]
  			  end

  			  action_node.elements.each("COMMENT") do |comm|
  				actions << ['comment', comm.attributes.find{|a| a[0].upcase == 'REF'}[1].strip()]
  			  end

  			  if(action_node.elements['SCORERANGE'] != nil and action_node.elements['SCORERANGE/USE_GLOBALRANGE'] != nil)
  				actions << ['scorerange','useglobalrange']
  			  elsif(action_node.elements['SCORERANGE'] != nil)
  				srange = action_node.elements['SCORERANGE'].text().strip().split(',').map{|e| e.gsub(/^[ \(\)\n]*/,'').gsub(/[ \(\)\n]*$/, '') }
  				srange.map! do |e|
  				  e =~ /\s*(\S+)\s*TO\s*(\S+)\s*=>\s*(\S+)\s*/
  				  [$1,$2,$3]
  				end
  				actions << ['scorerange', srange]
  			  end
  			end
  			rules << [condition, actions]
  		  end
  		  @drugs << [name, fullname, rules]
  		end

  		#MUTATION_COMMENTS
  		dom.elements.each("*/MUTATION_COMMENTS/GENE") do |gene_node|
  		  gene_name = gene_node.elements['NAME'].text().strip()
  		  rules = []

  		  gene_node.elements.each("RULE") do |rule_node|
  			condition = rule_node.elements['CONDITION'].text().strip().gsub(/\s+/, ' ')
  			actions = []

  			rule_node.elements.each("ACTIONS") do |action_node|
  			  action_node.elements.each("LEVEL") do |level|
            actions << ['level', level.text().to_i]
  			  end

  			  action_node.elements.each("COMMENT") do |comm|
  				actions << ['comment', comm.attributes.find{|a| a[0].upcase == 'REF'}[1].strip()]
  			  end

  			  if(action_node.elements['SCORERANGE'] != nil and action_node.elements['SCORERANGE/USE_GLOBALRANGE'] != nil)
  				actions << ['scorerange','useglobalrange']
  			  elsif(action_node.elements['SCORERANGE'] != nil)
  				srange = defs.elements['SCORERANGE'].text().strip().split(',').map{|e| e.gsub(/^[ \(\)\n]*/,'').gsub(/[ \(\)\n]*$/, '') }
  				srange.map! do |e|
  				  e =~ /\s*(\S+)\s*TO\s*(\S+)\s*=>\s*(\S+)\s*/
  				  [$1,$2,$3]
  				end
  				actions << ['scorerange', srange]
  			  end
  			end
  			rules << [condition, actions]
  		  end

  		  @mutation_comments << [gene_name, rules]
  		end

  	end
  end

  #BNF Parsing
  def interp_condition(cond, aaseq)
    bnf = bnf_statement(cond + '|', aaseq)
    if(bnf.cond == false)
      puts "---Could not parse algorithm condition:  " + cond
      exit()
      return nil
    elsif(bnf.truth)
      #print "RESULT IS TRUE, SCORE: " + str(bnf.score)
      return [true, bnf.score, bnf.flags.uniq()]
    else
      #print "RESULT IS FALSE, SCORE: " + str(bnf.score)
      return [false, bnf.score, bnf.flags.uniq()]
    end
  end

  # booleancondition | scorecondition
  def bnf_statement(cond, aaseq)
    puts "statement: " + cond  if(@debug)

    [method(:bnf_booleancondition), method(:bnf_scorecondition)].each do |func|
      bnf = func.call(cond, aaseq)
      return bnf if(bnf.cond)
    end
    return BNFVal.new(false)
  end

  # condition condition2*;
  def bnf_booleancondition(cond, aaseq)
    puts "booleancondition: " + cond if(@debug)

    bnflist = []
    bnf = bnf_condition(cond, aaseq)
    if(bnf.cond)
      bnflist << bnf
      while(true)
        bnf = bnf_condition2(bnf.cond, aaseq)
        if(bnf.cond)
          bnflist << bnf
        else
          break
        end
      end
      #Use the logic
      left_truth = nil
      bnflist.each do |bnf|
        #print "boolean_truth " + str(bnf.truth)
        if(left_truth == nil)
          left_truth = bnf.truth
        else #Not quite proper logic order, but I don't think anybody is randomly mixing OR's and AND's so it should be okay
          if(bnf.logic == 'AND' and left_truth and bnf.truth)
            left_truth = true
          elsif(bnf.logic == 'AND')
            left_truth = false
          elsif(bnf.logic == 'OR' and (left_truth or bnf.truth))
            left_truth = true
          elsif(bnf.logic == 'OR')
            left_truth = false
          end
        end
      end
      return BNFVal.new(bnflist[-1].cond, left_truth)
    else
      return BNFVal.new(false)
    end
  end

  #l_par booleancondition r_par | residue | excludestatement | selectstatement
  def bnf_condition(cond, aaseq)
    puts "condition: " + cond  if(@debug)
    [method(:bnf_booleancondition), method(:bnf_residue), method(:bnf_excludestatement), method(:bnf_selectstatement)].each do |func|
      if(func == method(:bnf_booleancondition))
        lpi = -1
        rpi = -1
        cnt = 0
        #immediaate break if it doesn't match this regexp:
#        next if (not re.match('^\s*\(', cond)) #TODO TODO TODO
        next if(not cond =~ /^\s*\(/)

        #Search for parens
        0.upto(cond.size() - 1) do |i|
          if(cond[i,1] == '(')
            lpi = i if(lpi == -1)
            cnt += 1
          elsif(cond[i,1] == ')')
            cnt -= 1
            if(cnt == 0)
              rpi = i
              break
            end
          else
            next
          end
        end

        next if(lpi == -1 or rpi == -1)

        bnf = func.call(cond[lpi + 1 ... rpi] + ' |', aaseq)
        return bnf if(bnf.cond)
      else
        bnf = func.call(cond, aaseq)
        return bnf if(bnf.cond)
      end
    end

    return BNFVal.new(false)
  end

  #logicsymbol condition;
  def bnf_condition2(cond, aaseq)
    puts "condition2: " + cond if(@debug)
    bnf_logic = bnf_logicsymbol(cond, aaseq)
    if(bnf_logic.cond)
      bnf = bnf_condition(bnf_logic.cond, aaseq)
      if(bnf.cond)
        bnf.logic = bnf_logic.logic
        return bnf
      end
    end
    return BNFVal.new(false)
  end

  #and | or
  def bnf_logicsymbol(cond, aaseq)
    puts "logicsymbol: " + cond if(@debug)
    if(cond =~ /^\s*AND\s*/i)
      bnf = BNFVal.new(cond.sub(/^\s*AND\s*/, ''))
      bnf.logic = 'AND'
      return bnf
    elsif(cond =~ /^\s*OR\s*/i)
      bnf = BNFVal.new(cond.sub(/^\s*OR\s*/, ''))
      bnf.logic = 'OR' #I think this works????
      return bnf
    end
    return BNFVal.new(false)
  end

  #[originalaminoacid]:amino_acid? integer [mutatedaminoacid]:amino_acid+ |
  #not [originalaminoacid]:amino_acid? Integer [mutatedaminoacid]:amino_acid+ |
  #[originalaminoacid]:amino_acid? integer l_par not [mutatedaminoacid]:amino_acid+ r_par
  def bnf_residue(cond, aaseq) #this looks hard yo
    puts "residue: " + cond if(@debug)
    truth = false
    #I think we'll have to go the regexp route here.  Haha.
    #mo_a = re.search('^\s*[ARNDCEQGHILKMFPSTWYVid]?\s*(\d+)\s*([ARNDCEQGHILKMFPSTWYVid]+)', cond)
    #mo_b = re.search('^\s*NOT\s*[ARNDCEQGHILKMFPSTWYVid]?\s*(\d+)\s*([ARNDCEQGHILKMFPSTWYVid]+)', cond)
    #mo_c = re.search('^\s*[ARNDCEQGHILKMFPSTWYVid]?\s*(\d+)\s*\(\s*NOT\s*([ARNDCEQGHILKMFPSTWYVid]+)\s*\)', cond)


	if(cond =~ /^\s*TRUE\s*/) #new special case
	  puts 'mo_true' if(@debug)
      bnf = BNFVal.new(cond.sub(/^\s*TRUE\s*/, ''), true)
      puts "Residue: " + bnf.truth.to_s if(@debug)
      return bnf
    elsif(cond =~ /^\s*[ARNDCEQGHILKMFPSTWYVid]?\s*(\d+)\s*([ARNDCEQGHILKMFPSTWYVid]+)/)
      puts 'mo_a' if(@debug)
      loc = $1.to_i
      aas = $2
      if(aaseq.size() > loc)
        aaseq[loc - 1].each do |aa|
          truth = true if(aas.include?(aa))
        end
      end
      bnf = BNFVal.new(cond.sub(/^\s*[ARNDCEQGHILKMFPSTWYVid]?\s*(\d+)\s*([ARNDCEQGHILKMFPSTWYVid]+)/, ''), truth)
      puts "Residue: " + bnf.truth.to_s if(@debug)
      return bnf
    elsif(cond =~ /^\s*NOT\s*[ARNDCEQGHILKMFPSTWYVid]?\s*(\d+)\s*([ARNDCEQGHILKMFPSTWYVid]+)/)
      puts 'mo_b' if(@debug)
      loc = $1.to_i
      aas = $2
      truth = true
      if(aaseq.size() > loc)
        aaseq[loc - 1].each do |aa|
          truth = false if(aas.include?(aa))
        end
      else
        truth = false # ????UNKNOWN
      end
      bnf = BNFVal.new(cond.sub(/^\s*NOT\s*[ARNDCEQGHILKMFPSTWYVid]?\s*(\d+)\s*([ARNDCEQGHILKMFPSTWYVid]+)/, ''), truth)
      puts "Residue: " + bnf.truth.to_s if(@debug)
      return bnf
#=begin
	elsif(cond =~ /^\s*[ARNDCEQGHILKMFPSTWYVid]?\s*(\d+)\s*\!\s*([ARNDCEQGHILKMFPSTWYVid]+)/)
	  puts 'mo_b(2)' if(@debug)  #Slightly different than NOT, more of a NONE.
      loc = $1.to_i
      aas = $2
      truth = true
      if(aaseq.size() > loc)
        truth = false if(aaseq[loc - 1].all?() {|aa| aas.include?(aa) })
      else
        truth = false # ????UNKNOWN
      end
      bnf = BNFVal.new(cond.sub(/^\s*[ARNDCEQGHILKMFPSTWYVid]?\s*(\d+)\s*\!\s*([ARNDCEQGHILKMFPSTWYVid]+)/, ''), truth)
      puts "Residue: " + bnf.truth.to_s if(@debug)
      return bnf
#=end
    elsif(cond =~ /^\s*[ARNDCEQGHILKMFPSTWYVid]?\s*(\d+)\s*\(\s*NOT\s*([ARNDCEQGHILKMFPSTWYVid]+)\s*\)/)
      puts 'mo_c' if(@debug)
      loc = $1.to_i
      aas = $2
      truth = true
      if(aaseq.size() > loc)
        aaseq[loc - 1].each do |aa|
          if(aas.include?(aa))
            truth = false
          elsif(aa != '*')
            truth = true
            break
          end
        end
      else
        truth = false # ????UNKNOWN
      end
      bnf = BNFVal.new(cond.sub(/^\s*[ARNDCEQGHILKMFPSTWYVid]?\s*(\d+)\s*\(\s*NOT\s*([ARNDCEQGHILKMFPSTWYVid]+)\s*\)/, ''), truth)
      puts "Residue: " + bnf.truth.to_s if(@debug)
      return bnf
    end
    return BNFVal.new(false)
  end

  #exclude residue
  def bnf_excludestatement(cond, aaseq)
    puts "excludestatement: " + cond if(@debug)
    if(cond =~ /^\s*EXCLUDE\s*/i)
      bnf = bnf_residue(cond.sub(/^\s*EXCLUDE\s*/, ''), aaseq)
      if(bnf.cond)
        bnf.truth = !bnf.truth
        return bnf
      end
    end
    return BNFVal.new(false)
  end

  #select selectstatement2
  def bnf_selectstatement(cond, aaseq)
    puts "selectstatement: " + cond if(@debug)
    if(cond =~ /^\s*SELECT\s*/i)
      bnf = bnf_selectstatement2(cond.sub(/^\s*SELECT\s*/, ''), aaseq)
      if(bnf.cond)
        return bnf
      end
    end
    return BNFVal.new(false)
  end

  #exactly integer from l_par selectlist r_par |
  #atleast integer from l_par selectlist r_par |
  #notmorethan integer from l_par selectlist r_par |
  #atleast [atleastnumber]:integer logicsymbol notmorethan [notmorethannumber]:integer from l_par selectlist r_par
  def bnf_selectstatement2(cond, aaseq)
    puts "selectstatement2: " + cond if(@debug)
    lparen = -1
    rparen = -1
    cnt = 0
    0.upto(cond.size() - 1) do |i|
      if(cond[i,1] == '(' and cnt == 0)
        lparen = i
        cnt += 1
      elsif(cond[i,1] == '(')
        cnt += 1
      elsif(cond[i,1] == ')' and cnt == 1)
        rparen = i
        cnt -= 1
        break
      elsif(cond[i,1] == ')')
        cnt -= 1
      end
    end

    if(lparen == -1 or rparen == -1)
      return BNFVal.new(false)
    end

    #get the list items
    bnflist = bnf_selectlist(cond[lparen+1 ... rparen] + ' |', aaseq)

    cond =~ /^\s*(EXACTLY\s*(\d+)|ATLEAST\s*(\d+)\s*(AND|OR)\s*NOTMORETHAN\s*(\d+)|ATLEAST\s*(\d+)|NOTMORETHAN\s*(\d+))\s*from\s*\(\s*(.+)\s*\)/i
    mo_a = [$0, $1, $2, $3, $4, $5, $6, $7]
    #mo_a = re.search('^\s*(EXACTLY\s*(\d+)|ATLEAST\s*(\d+)\s*(AND|OR)\s*NOTMORETHAN\s*(\d+)|ATLEAST\s*(\d+)|NOTMORETHAN\s*(\d+))\s*from\s*\(\s*(.+)\s*\)', cond, flags=re.I)

    atleastn = -1
    atmostn = -1
    exactlyn = -1
    cnt = 0
    logic = nil
    type = mo_a[1] #mo_a.group(1)

    bnflist.each do |bnf|
      cnt += 1 if(bnf.cond and bnf.truth)
    end

    puts "cnt: #{cnt}" if(@debug)
    #if re.search('^\s*ATLEAST\s*(\d+)\s*(AND|OR)\s*NOTMORETHAN', type, flags=re.I):
    if(type =~ /^\s*ATLEAST\s*(\d+)\s*(AND|OR)\s*NOTMORETHAN/i)
      atleastn = mo_a[3].to_i
      logic = mo_a[4]
      atmostn = mo_a[5].to_i
      puts "atleast: " + atleastn.to_s + ", atmost: " + atmostn.to_s if(@debug)
      if(cnt >= atleastn and cnt <= atmostn)
        return BNFVal.new(cond[rparen + 1 .. -1], true)
      else
        return BNFVal.new(cond[rparen + 1 .. -1], false)
      end
    elsif(type =~ /^\s*EXACTLY/i)
      exactlyn = mo_a[2].to_i
      puts "exactly: " + exactlyn.to_s if(@debug)
      if(cnt == exactlyn)
        return BNFVal.new(cond[rparen + 1 .. -1], true)
      else
        return BNFVal.new(cond[rparen + 1 .. -1], false)
      end
    elsif(type =~ /^\s*ATLEAST/i)
      atleastn = mo_a[6].to_i
      puts "atleast: " + atleastn.to_s  if(@debug)
      if(cnt >= atleastn)
        return BNFVal.new(cond[rparen + 1 .. -1], true)
      else
        return BNFVal.new(cond[rparen + 1 .. -1], false)
      end
    elsif(type =~ /^\s*NOTMORETHAN/i)
      atmostn = mo_a[7].to_i
      puts "atmost: " + atmostn.to_s if(@debug)
      if(cnt <= atmostn)
        return BNFVal.new(cond[rparen + 1 .. -1], true)
      else
        return BNFVal.new(cond[rparen + 1 .. -1], false)
      end
    end

    return BNFVal.new(false)
  end


  #residue listitems*
  def bnf_selectlist(cond, aaseq)
    puts "selectlist: " + cond if(@debug)
    bnflist = []
    bnf = bnf_residue(cond, aaseq)
    if(bnf.cond)
      newcond = bnf.cond
      bnflist << bnf
      while(true)
#        tmp = re.search('^\s*,\s*', newcond) #doesn't seem to be used.
        newcond = newcond.sub(/^\s*,\s*/, '')
        bnf = bnf_residue(newcond, aaseq)
        if(bnf.cond)
          bnflist << bnf
          newcond = bnf.cond
        else
          break
        end
      end
      return bnflist
    end
    return [BNFVal.new(false)]
  end

  #score from l_par scorelist r_par
  def bnf_scorecondition(cond, aaseq)
    puts "bnf_scorecondition: " + cond if(@debug)

    if(cond =~ /^\s*score\s*from\s*\((.+)\)\s*|/i and $1)
      score = 0.0
	  flags = []
      bnf_list = bnf_scorelist($1 + '|', aaseq)
      if(bnf_list[0].cond)
        newcond = cond
        bnf_list.each do |bnf|
          if(bnf.cond)
            newcond = bnf.cond
            score += bnf.score
			flags += bnf.flags
          else
            break
          end
        end
      else
        return BNFVal.new(false)
      end
      #I guess evalate the truth values?
      return BNFVal.new(newcond, false, score=score,flags=flags.uniq())
      #return bnf
    end
    return BNFVal.new(false)
  end


  #scoreitem scoreitems*
  def bnf_scorelist(cond, aaseq)
    puts "bnf_scorelist: " + cond  if(@debug)
    bnflist = []
    bnf = bnf_scoreitem(cond, aaseq)
    puts "S:" + bnf.score.to_s if(@debug)
    if(bnf.cond) #***
      bnflist << bnf
      newcond = bnf.cond
      while(newcond)
        #check for comma
        if(newcond =~ /^\s*,\s*/)
          newcond = newcond.sub(/^\s*,\s*/, '')
          bnf = bnf_scoreitem(newcond, aaseq)
          break if(bnf.cond == false)
          puts "S:" + bnf.score.to_s if(@debug)
          newcond = bnf.cond
          bnflist << bnf
        else
          break
        end
      end
      return bnflist
    end
    return [BNFVal.new(false)]
  end

  #booleancondition mapper min? number |
  #max l_par scorelist r_par
  def bnf_scoreitem(cond, aaseq)
    puts "bnf_scoreitem: " + cond  if(@debug)
    #Trickier than we think, as the subexpressions could have parens.  Need to do the counting game.  Sadly.
    if(cond =~ /^\s*MAX\s*\(/i)
      #print "to mob, or not to mob?"
      lparen = -1
      rparen = -1
      cnt = 0
      0.upto(cond.size() - 1) do |i|
        if(cond[i,1] == '(' and cnt == 0)
          lparen = i
          cnt += 1
        elsif(cond[i,1] == '(')
          cnt += 1
        elsif(cond[i,1] == ')' and cnt == 1)
          rparen = i
          cnt -= 1
          break
        elsif(cond[i,1] == ')')
          cnt -= 1
        end
      end

      if(lparen == -1 or rparen == -1)
        return BNFVal.new(false)
      end
      newcond = cond[lparen + 1 ... rparen]
      bnflist = bnf_scorelist(newcond, aaseq)
      score = -999 #close enough to infinity.

      if(bnflist[0].cond)
        bnflist.each do |bnf|
          if(bnf.cond)
            newcond = bnf.cond
            if(bnf.score > score and bnf.score != 0.0)
              score = bnf.score
            end
          end
        end
        score = 0.0 if(score == -999)
        return BNFVal.new(cond[rparen+1 .. -1], false, score=score)
      end
    elsif(cond =~ /^\s*([^=>]+)\s*=>\s*(min)?\s*(-?\d+\.?\d*)\s*(.+)$/i)
      #mo_a has 4 groups, the booleanconditon, an optional pointless MIN, and the score, and then the rest of the string
      #I think we need to match a dash for negative numbers yo
      mo_a = [$0, $1, $2, $3, $4]
      bnf = bnf_booleancondition(mo_a[1], aaseq)
      bnf_score = mo_a[3].to_f
      if(bnf.cond)
        if(bnf.truth)
          bnf.score = bnf_score
        end
        bnf.cond = mo_a[4]
        return bnf
	  end
	elsif(cond =~ /^\s*([^=>]+)\s*=>\s*"([^"]+)"\s*(.+)$/i) #Direct interp call, adds flag.
      #mo_a has 3 groups, the booleanconditon, and the interp, and then the rest of the string
      #I think we need to match a dash for negative numbers yo
      mo_a = [$0, $1, $2, $3]
      bnf = bnf_booleancondition(mo_a[1], aaseq)
      bnf_score = 0
	  bnf_flag = mo_a[2]
      if(bnf.cond)
        if(bnf.truth)
          bnf.score = bnf_score
		  bnf.flags << bnf_flag
        end
        bnf.cond = mo_a[3]
        return bnf
	  end
    end
    return BNFVal.new(false)
  end

  #handles a couple undocumented comment filtering things.
  def comment_filter(comment, aaseq, region=nil)
    #listMutsIn
    #tmp = re.match('^.*\$listMutsIn\{([^\}]+)\}.*$', comment)
    #comment =~ /^.*\$listMutsIn\{([^\}]+)\}.*$/
    #tmp = $1
    if(comment =~ /^.*\$listMutsIn\{([^\}]+)\}.*$/)
      tmp = $1
      tmporig = tmp.to_s
      tmporig = tmporig.gsub('(', '\(').gsub(')','\)')
      muts = tmp.split(',')
      final = []
      muts.each do |mut|
        match = ''
        loc = ''
        #If it matches \d+\(NOT \w+\), then we got to do something fancy
        if(mut =~ /[a-z]?(\d+)\(NOT\s+([a-z]+)\)/i)
          tmpa = [$0, $1, $2]
          loc = tmpa[1]
          match = 'ARNDCEQGHILKMFPSTWYVid'
          tmpa[2].split('').each do |ch|
            match = match.gsub(ch, '')
          end
        elsif(mut =~ /[a-z]?(\d+)([a-z]+)/i)
          tmpb = [$0, $1, $2]
          loc = tmpb[1]
          match = tmpb[2]
        end
        aas = aaseq[loc.to_i - 1]
        subs = ''
        aas.each do |aa|
          if(match.include?(aa))
            subs += aa
          end
        end
        if(subs != '')
          subs = subs.split('')
          subs.sort!()
          subs = subs.join('')
          if(region == 'PR')
            final << (@pr_std[loc.to_i - 1, 1] + loc + subs)
          elsif(region == 'RT')
            final << (@rt_std[loc.to_i - 1, 1] + loc + subs)
          elsif(region == 'IN')
            final << (@int_std[loc.to_i - 1, 1] + loc + subs)
          else
            final << (loc + subs)
          end
        end
      end

      comment = comment.sub(/\$listMutsIn\{#{tmporig}\}/, final.join(', '))
      comment = comment.sub(/ \(\)/, '') #get rid of empty brackets.
    end

    #numberOfMutsIn
    if(comment =~ /^.+\$numberOfMutsIn\{([^\}]+)\}.+$/)
      tmpmatch = $1
      muts = tmpmatch.split(',')
      cnt = 0
      muts.each do |mut|
        mut =~ /^(\d+)([a-z]+)$/i
        tmp = [$0, $1, $2]
        aas = aaseq[tmp[1].to_i - 1]
        aas.each do |aa|
          if(tmp[2].include?(aa))
            cnt += 1
          end
        end
      end
      comment = comment.sub(/\$numberOfMutsIn\{#{tmpmatch}\}/, cnt.to_s)
    end

    comment = comment.sub(/  /, ' ') #Make spacing more like sierra
    #unicod = '\uf0b1'  #blorp, how to do this in ruby?
    #unicode = ['F','0','B','1'].pack('U')
    #Maybe ignore the unicode?????
    #comment = comment.sub(unicod, '+/-') #Fixing crazy unicode characters
    #TODO:  If things break, you might have to fix stuff here.
    return comment
  end



  #Most important method.
  def interpret(aaseq, region)
    result = AsiResult.new()
    result.alg_name = @alg_name
    result.alg_version = @alg_version
    genes = @gene_def.find_all{|e| e[0] == region}[0][1]
    drs = @drug_class.find_all{|e| genes.include?(e[0]) }
    #print drs
    drs.each do |drcls|
      cls = drcls[0]
      drcls[1].each do |drname|
        drug = @drugs.find_all {|e| e[0] == drname}[0]
        drug_result = AsiDrugResult.new()
        drug_result.code = drug[0]
        drug_result.name = drug[1]
        drug_result.drug_class = cls

        drug[2].each do |rule|
          cond = rule[0]
          actions = rule[1]
          interp = interp_condition(cond, aaseq)

          if(interp)
            score = interp[1]
            truth = interp[0]

            puts "Final Score:" + score.to_s if(@debug)
            puts "Final Score:" + truth.to_s  if(@debug)
            next if(truth == false and score == 0.0)

            actions.each do |act|
              if(act[0] == 'level')
                if(act[1].to_i > drug_result.level)
                  drug_result.level = act[1].to_i
                end
              elsif(act[0] == 'comment')
                comm = act[1] #why is this even here.
                comm =@comment_def.find_all{|e| e[0] == act[1]}[0]
                comment = comm[1]
                #while(re.search('\$numberOfMutsIn{', comment) or re.search('\$listMutsIn{', comment))
                while(comment =~ /\$numberOfMutsIn\{/ or comment =~ /\$listMutsIn\{/)
                  comment = comment_filter(comment, aaseq, region)
                end
                drug_result.comments << comment
              elsif(act[0] == 'scorerange')
                drug_result.score = score
                scorerange = act[1]
                if(scorerange == 'useglobalrange')
                  scorerange = @global_range
                end

                scorerange.each do |rng|
                  if(rng[0] == '-INF')
                    rng[0] = -99999  #that is close enough to negative infinity.
                  else
                    begin
                      rng[0] = rng[0].to_f
                    rescue  #Wat?
                      rng[0] = rng[0].to_f
                    end
                  end
                  if(rng[1] == 'INF')
                    rng[1] = 99999 #that is close enough to infinity.
                  else
                    begin
                      rng[1] = rng[1].to_f
                    rescue #WAT?
                      rng[1] = rng[1].to_f
                    end
                  end

                  if(drug_result.score >= rng[0] and drug_result.score <= rng[1])
                    if(rng[2].to_i > drug_result.level)
                      drug_result.level = rng[2].to_i
                    end
                    break
                  end
                end
                if(drug_result.level == nil)
                  raise "drug score range level parsing error"
                end
              end
            end
          elsif(interp == nil)
            puts "ERROR in condition: " + cond
          end
        end
        result.drugs << drug_result
      end
    end

    result.drugs.sort!(){|a,b| a.drug_class == b.drug_class ?  a.code <=> b.code : b.drug_class <=> a.drug_class}
    #result.drugs.sort(key=lambda e: e.code)
    #result.drugs.sort(key=lambda e: e.drug_class, reverse=true)

    #comments
    @mutation_comments.each do |gene|
      next if(gene[0] != region)

      gene[1].each do |mut|
        cond = mut[0]
        actions = mut[1]

        interp = interp_condition(cond, aaseq)
        if(interp and interp[0])
          actions.each do |act|
            comm = @comment_def.find_all(){|e| e[0] == act[1]}[0]
            comment = comment_filter(comm[1], aaseq, region)
            mut = comm[0]
            #puts mut
            mut =~ /#{region}(\d+)/
            result.mutations << [$1.to_i, aaseq[$1.to_i + 1]]
            result.mutation_comments << comment
          end
        elsif(interp == nil)
          puts "ERROR in condition: " + cond
        end
      end
    end

    #mutations
    #~ for mut in self.mutations:
      #~ if mut[0] == region:
        #~ aas = aaseq[int(mut[1]) - 1]
        #~ subs = ''
        #~ for aa in aas:
          #~ if aa in mut[2]:
            #~ subs += aa
        #~ if subs != '':
          #~ subs = list(subs)
          #~ subs.sort()
          #~ subs = ''.join(subs)
          #~ moot = ''
          #~ if region == 'PR':
            #~ moot = self.pr_std[int(mut[1]) - 1] + mut[1] + subs
          #~ elif region == 'RT':
            #~ moot = self.rt_std[int(mut[1]) - 1] + mut[1] + subs
          #~ elif region == 'IN':
            #~ moot = self.int_std[int(mut[1]) - 1] + mut[1] + subs
          #~ else:
            #~ moot = loc + subs

          #~ if mut[3] == 'Other':
            #~ result.mutations_other.append(moot)
          #~ else:
            #~ result.mutations.append(moot)


    return result
  end
end


class AsiDrugResult
  attr_accessor :name, :code, :drug_class, :score, :level, :comments
  def initialize()
    @name = ''
    @code = ''
    @drug_class = ''
    @score = 0.0
    @level = 1
    @comments = []  #Do we need this?
  end

end

class AsiResult
  attr_accessor :alg_name, :alg_version, :mutation_comments, :mutations, :mutations_other, :drugs
  def initialize()
    @alg_name = ''
    @alg_version = ''

    @mutation_comments = []
    @mutations = []  #mut hash yo!  Don't forget yer proteins!  (Only used for HIVDB, and its a bit of a hack)
    @mutations_other = []  #mut hash yo!  Don't forget yer proteins!  (Only used for HIVDB, and its a bit of a hack)
    @drugs = []  #drug hash yo!  Don't forget your scores!
  end
end

=begin
#Testing
AA_HASH = {
      'ttt' => 'F', 'tct' => 'S', 'tat' => 'Y', 'tgt' => 'C',
      'ttc' => 'F', 'tcc' => 'S', 'tac' => 'Y', 'tgc' => 'C',
      'tta' => 'L', 'tca' => 'S', 'taa' => '*', 'tga' => '*',
      'ttg' => 'L', 'tcg' => 'S', 'tag' => '*', 'tgg' => 'W',

      'ctt' => 'L', 'cct' => 'P', 'cat' => 'H', 'cgt' => 'R',
      'ctc' => 'L', 'ccc' => 'P', 'cac' => 'H', 'cgc' => 'R',
      'cta' => 'L', 'cca' => 'P', 'caa' => 'Q', 'cga' => 'R',
      'ctg' => 'L', 'ccg' => 'P', 'cag' => 'Q', 'cgg' => 'R',

      'att' => 'I', 'act' => 'T', 'aat' => 'N', 'agt' => 'S',
      'atc' => 'I', 'acc' => 'T', 'aac' => 'N', 'agc' => 'S',
      'ata' => 'I', 'aca' => 'T', 'aaa' => 'K', 'aga' => 'R',
      'atg' => 'M', 'acg' => 'T', 'aag' => 'K', 'agg' => 'R',

      'gtt' => 'V', 'gct' => 'A', 'gat' => 'D', 'ggt' => 'G',
      'gtc' => 'V', 'gcc' => 'A', 'gac' => 'D', 'ggc' => 'G',
      'gta' => 'V', 'gca' => 'A', 'gaa' => 'E', 'gga' => 'G',
      'gtg' => 'V', 'gcg' => 'A', 'gag' => 'E', 'ggg' => 'G',
}

AMBIG = { 'A'=> ['A'], 'G'=> ['G'], 'T'=> ['T'], 'C'=> ['C'],
'R' => ['A', 'G'], 'Y' => ['C', 'T'], 'K' => ['G', 'T'],
'M' => ['A', 'C'], 'B' => ['C', 'G', 'T'], 'D' => ['A', 'G', 'T'],
'H' => ['A', 'C', 'T'], 'V' => ['A', 'C', 'G'], 'S' => ['C', 'G'],
'W' => ['A', 'T'], 'N' => ['A', 'C', 'G', 'T'], 'X' => ['X'], '-' => ['-'] }


#Generates all possible non-ambigious nucleotides from an ambiguous nucleotide
def generate(nuc)
	posa = AMBIG[nuc[0,1]]
	posb = AMBIG[nuc[1,1]]
	posc = AMBIG[nuc[2,1]]

	if(nuc == 'XXX')
		return ['X']
	end
	nuclist = []
	posa.each do |a|
		posb.each do |b|
			posc.each do |c|
				nuclist.push(a + b + c)
			end
		end
	end
	return nuclist
rescue StandardError => error
	puts error
	puts nuc
	return nil
end

#returns an array of arrays,
#translate_complete_to_array('AAARRR') => [["K"], ["K", "R", "E", "G"]]
def translate_complete_to_array(sequence)
	protseq = [];
	0.upto((sequence.size() / 3) - 1) do |i|
		list = generate(sequence[i * 3, 3])
		list.map! { |entry|
				if(entry == 'X')
					'X'
				else
					raw_translator(entry)
				end
		}.uniq!

		if(list.size > 1)
			tmp = []
			list.each do |entry|
				tmp.push(entry)
			end

			protseq.push(tmp)
		elsif(list.size == 1)
			protseq.push([list[0]])
		else
			protseq.push(['*'])
		end
	end

	return protseq;
end

def raw_translator(str, keepdashs = false)
	aa = ''
	0.upto(((str.length - 1) / 3)) do |i|
        if(keepdashs and str[i * 3, 3] == '---')
            aa += '-'
            next
        end
		x = AA_HASH[str[i * 3, 3].downcase]
		x = 'X' if(x == nil)
		aa += x
	end
	return aa;
end

int_seq = 'TTTTTAGATGGAATAGATAAGGCCCAAGATGAACATGAGAAATATCACAGTAATTGGAGAGCAATGGCTAGTGATTTTAACCTGCCACCTGTAGTAGCAAAAGAAATAGTAGCCAGCTGTGATAAATGTCAGCTAAAAGGAGAAGCCATGCATGGACAAGTAGACTGTAGTCCAGGAATATGGCAACTAGATTGTACACATTTAGAAGGAAAAGTTATCATGGTAGCAGTTCATGTAGCCAGTGGATATATAGAAGCAGAAGTTATTCCAGCAGAAACAGGGCAGGAAACAGCATATTTTCTTTTAAAATTAGCAGGAAGATGGCCGGTAAAAACAATACATACAGACAATGGCAGCAATTACACCAGTGCTACGGTTAAGGCCGCCTGTTGGTGGGCGGGAATCAAGCAGGAATTTGGAATTCCCTACAATCCCCAAAGTCAAGGAGTAATAGAATCTATGAATAAAGAATTAAAGAAAATTATAAGACAGGTAAGAGATCAGGCTGAACATCTTAAGACAGCAGTACAAATGGCAGTATTCATCCACAATTTTAAAAGAAAAGGGGGGATTGGGGGGTACAGTGCAGGGGAAAGAATAGTAGACATAATAGCAACAGACATACAAACTAAAGAATTACAAAAACAAATTACAAAAATTCAAAATTTTCGGGTTTATTACAGGGACAGCAGAAATCCACTTTGGAAAGGACCAGCAAAGCTCCTCTGGAAAGGTGAAGGGGCAGTAGTAATACAAGATAATAGTRACATAAAAGTAGTGCCAAGAAGAAAAGCAAAGATCATTAGGGATTATGGAAAACAGATGGCAGGTGATGATTGTGTGGCAAGTAGACAGGATGAGGAT'
pr_seq = 'CCTCAAATCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGRGGGCAGCTAAMGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGACATGGARTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGRTACCCATAGAAATTTGTGGACAYAAAACTATAGGTWCAGTATTAATAGGACCTACACCWGTTAACATAATTGGAAGAAATCTGATGAYTCAGCTTGGTTGCACTTTAAATTTT'
rt_seq = 'CCCATTAGTCCTATTGAAACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAGGTYAARCAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAAATGGAAAAGGAAGGRAAGATTTCAAAAATTGGACCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAARAGAACTCAAGACTTCTGGGAAGTTCAATTAGGAATACCACATCCYGCAGGGTTAAAAAAGAAMAAGTCAGTAACAGTACTRGATGTGGGTGATGCATATTTTTCAGTTCCCTTATATGAAGACTTCAGGAAGTATACTGCATTCACCATACCTAGYACAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTGCCACAAGGATGGAAAGGATCACCAGCAATATTCCAAAGTAGCATGATAAAAATCTTAGAGCCTTTCAGAAAACAAAATCCAGARATAGTCATCTATCAATACGTGGATGATTTGTATGTAGSATCTGACTTAGAAATAGGGCAGCATAGAACAAAGATAGAGGAACTGAGAGCACATCTRTTRAAGTGGGGATTTACCACACCAGACAAAAAACATCAGAAAGAGCCTCCATTCCTTTGGATGGGTTATGAACTCCATCCTGATAAATGGACR'

aa_seq = translate_complete_to_array(rt_seq)


algorithm_hivdb = AsiAlgorithm.new(nil)
#algorithm_hivdb.debug = true
#algorithm_anrs = AsiAlgorithm.new("ANRS_max.xml")
#algorithm_rega = AsiAlgorithm.new("RegaInst_max.xml")


res = algorithm_hivdb.interpret(aa_seq, 'RT')

res.drugs.each do |drug|
  puts drug.code + ":  SCORE: " + drug.score.to_s + " , LEVEL: " + drug.level.to_s
  drug.comments.each do |com|
    puts com
  end
end
puts "Mutation Comments:"
res.mutation_comments.each do |com|
  puts com
end


=end
