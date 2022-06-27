require 'alignment'

#tmp = align_it("ATGATAAATATTATACGCGTCA", "ATTATA", 3, 1)
tmp = align_it("ATGATAAATATTATACGC$$$A", "ATTATATAG", 3, 1)
puts tmp[0]
puts tmp[1]
puts "WOO!"
gets