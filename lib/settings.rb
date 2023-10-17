#interprets the config/settings.txt file.


class Settings
  @@settings = {}

  @@definitions = {
    'debug' => {type: 'Boolean', default: false, }, #test this
    'debug-genes' => {type: 'String', default: "pr,rt,int", },
    'debug-seq-limit' => {type: 'Integer', default: 999999, },

    'output-filter-direction' => {type: 'Boolean', default: false,},

    'jobfile-update-interval' => {type: 'Integer', default: 15, },

    'mix-threshold' => {type: 'Float', default: 0.25, },
    'align-gap-init' => {type: 'Float', default: 0.25, },
    'align-gap-extend' => {type: 'Float', default: 0.25, },
    'reference-match-threshold' => {type: 'Float', default: 0.70, },

    'trim-gene-edges' => {type: 'Boolean', default: true, },
    'trim-edges-dash-threshold' => {type: 'Float', default: 0.20, },

    'homopolymer-fixes' => {type: 'Boolean', default: false, },
    'homopolymer-fixes-fill-gaps' => {type: 'Boolean', default: false, },
    'reduce-dual-homopolymers' => {type: 'Boolean', default: false}, #experimental

    'alignment-optimize' => {type: 'Boolean', default: true,},
    'alignment-optimize-pad-size' => {type: 'Integer', default: 6,},

    'reject-poor-aligned-ends' => {type: 'Boolean', default: false,},
    'reject-poor-aligned-ends-window' => {type: 'Integer', default: 15,},
    'reject-poor-aligned-ends-size-threshold' => {type: 'Float', default: 0.20,},
    'reject-poor-aligned-ends-match-threshold' => {type: 'Float', default: 0.65,},

    'sweep-inserts' => {type: 'Boolean', default: true, },

    'optimization-stop-poor-alignments' => {type: 'Boolean', default: true, },
    'optimization-stop-poor-alignments-threshold' => {type: 'Float', default: 0.10, },
    'optimization-clipped-alignment' => {type: 'Boolean', default: true, },
    'optimization-clipped-alignment-buffer' => {type: 'Integer', default: 60, },
    'optimization-coverage-limit' => {type: 'Boolean', default: true, },
    'optimization-coverage-limit-target' => {type: 'Integer', default: 1000, },
    'optimization-mixturize-subtypes' => {type: 'Boolean', default: true, },

    'region-def' => {type: 'String', default: "pr:1-99 rt:1-440 int:1-288", }, #new, TODO
  }

  def Settings.settings()
    return @@settings
  end

  def Settings.[](key)
    return @@settings[key.downcase()]
  end

  def Settings.load(filename)
    File.open(filename, 'r') do |file|
      file.each_line do |line|
        next if(line.strip() == '' or line[0] == '#')
        tmp = line.strip().split('=')

        #check if this is a valid setting.
        definition = @@definitions[tmp[0].downcase()]
        if(definition.nil?)
          #raise "Invalid setting:  #{tmp[0]} - Unknown setting."
        elsif(definition[:type] == "Integer" and Integer(tmp[1]).nil?)
          raise "Invalid setting:  #{tmp[0]} - Not an integer."
        elsif(definition[:type] == "Integer")
          @@settings[tmp[0].downcase()] = Integer(tmp[1])
        elsif(definition[:type] == "Float" and Float(tmp[1]).nil?)
          raise "Invalid setting:  #{tmp[0]} - Not an float."
        elsif(definition[:type] == "Float")
          @@settings[tmp[0].downcase()] = Float(tmp[1])
        elsif(definition[:type] == 'Boolean' and tmp[1].downcase() == 'true')
          @@settings[tmp[0].downcase()] = true
        elsif(definition[:type] == 'Boolean' and tmp[1].downcase() == 'false')
          @@settings[tmp[0].downcase()] = false
        elsif(definition[:type] == 'Boolean')
          raise "Invalid setting:  #{tmp[0]} - Must be true or false."
        elsif(definition[:type] == 'String')
          @@settings[tmp[0].downcase()] = tmp[1]
        else
          raise "Unknown setting type - #{definition[:type]}"
        end
      end
    end

    #add defaults.
    @@definitions.each do |key, val|
      if(@@settings[key].nil?)
        @@settings[key] = val[:default]
      end
    end

  end

end
