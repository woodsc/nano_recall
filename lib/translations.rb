#basic useful translation methods

AA_TRANS = {
  'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
  'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
  'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => '*',
  'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',

  'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
  'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
  'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
  'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',

  'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
  'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
  'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
  'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',

  'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
  'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
  'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
  'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G',
  '---' => '-',
}

#puts AA_TRANS.invert.keys.sort().map(){|a| "'#{a}'"}.join(',')

NUC_MIX = {
  'A'=> ['A'],
  'G'=> ['G'],
  'T'=> ['T'],
  'C'=> ['C'],
  'R' => ['A', 'G'].sort,
  'Y' => ['C', 'T'].sort,
  'K' => ['G', 'T'].sort,
  'M' => ['A', 'C'].sort,
  'S' => ['G', 'C'].sort,
  'W' => ['A', 'T'].sort,
  'B' => ['C', 'G', 'T'].sort,
  'D' => ['A', 'G', 'T'].sort,
  'H' => ['A', 'C', 'T'].sort,
  'V' => ['A', 'C', 'G'].sort,
  'N' => ['A', 'C', 'T', 'G'].sort,
  'X' => ['X']
}
NUC_MIX.default = ['X']

def nuc_to_aa(nuc)
  return AA_TRANS[nuc]
end

def nuc_to_aa_string(nuc)
  aa_str = ''
  0.upto((nuc.size() / 3) - 1) do |i|
    aa = AA_TRANS[nuc[i * 3, 3]]
    aa_str += aa.nil? ? '?' : aa
  end
  return aa_str
end

def nuc_to_aa_array(nuc)
  aa_str = []
  0.upto((nuc.size() / 3) - 1) do |i|
    aa_str << AA_TRANS[nuc[i * 3, 3]]
  end
  return aa_str
end

#turns quality values into a char string.
def qual_to_string(qual_array)
  str = ''
  qual_array.each do |qual|
    if(qual == nil or qual == '-')
      qualstr += '-'
    elsif((qual.to_i / 5).floor() < 50)
      qualstr += (qual.to_i / 5).floor().to_s
    else
      qualstr += '9'
    end
  end
  return str
end
