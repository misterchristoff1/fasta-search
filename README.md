# fasta-search

Usage: $ python3 fasta-search.py <mode> <input_file> ...

is - INTERRUPTION SEARCH
     Searches an input fasta file for possible one, two and three base 
     interruptions and gives as outputs: search sequence (optional); total 
     appearances of each interruption; number of reads in which each 
     interruption appears; allele frequency.
     Bases immediately preceding and following possible interruptions must be
     given as optional arguments, however these can be left blank by inputting
     empty quotes (''). 
     The preceding and following bases can be modified to search for 
     interruptions of greater than three bases, in this instance -v should 
     always be used for clarity. 

Usage: $ python3 fasta-search.py is <input_file> <preceding_bases> 
	        <following_bases> 
            
            optional: <minimum counts to display (default=1)> 
                      <[-rev] also search for reverse complements>
                      <[-v] display full search sequences>
                      <[-ord] sort the output into descending order by reads>
                      <[-o outfilename] outputs to file, if no file name is
                      given the default file name is [input_file].csv>
                        
ss - SEQUENCE SEARCH
     Searches an input fasta file for a given sequence and outputs: total reads
     in file; search sequence; total appearances of sequence; number of reads
     in which the sequence appears; allele frequency. 

     Usage: $ python3 fasta-search.py ss <input_file> <search_sequence>
         
            optional: <[-rev] also search for reverse complement>

sf - SEQUENCE FILTER
     Searches an input fasta file for a given sequence and outputs any reads 
     containing the sequence. 
     Use '> output.fasta' on the command line to write to a new file. 

     Usage: $ python3 fasta-search.py ss <input_file> <search_sequence>
            
              optional: <[-rev] also search for reverse complements>
