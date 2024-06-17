 
# TODO
[x] implement Translator class
    [] ensure translator is faster than BioSeq 
[] flesh out skeleton logic
[] implement orf finder (without orfanage)
    [] create test dataset with CHESS dataset
    [] translate in 6 frames and store all orfs above threshold
    [] multithread ? 
[] integrate psauron for scoring orfs
[] integrate mmseqs2 for clustering orfs
[] implement debug in verbose mode in stderr

# Notes
- seems like pip and conda both work, package in pip, development in conda?
    - will expect the user to install orfanage, mmseqs2 themself

- orf finding steps:
    - extract long open reading frames
    - identify orfs with homology to known proteins (optional)
    - predict likely coding regions based on homology
    - 


# Questions
- TransDecoder
    - is the output just one ORF, or like one for each transcript, or all ORFs, with a threshold?
    - how does td perform ORF search? 
        - what is the 6 frame finding? 

- mmseqs2
    - what is the utility of using this to cluster?
    - what will the memory overhead be / is there a way to reduce for purpose of this tool? 

- development
    - what is the best method for packaging and shipping this pipeline?

- orf finder
    - is it worth it to encode the dna sequence into 2 bit to make it marginally faster? -> benchmark

# Links
https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=tgencodes#SG33
https://github.com/TransDecoder/TransDecoder/wiki