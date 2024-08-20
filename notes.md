 
# TODO
[x] implement Translator class
    [] ensure translator is faster than BioSeq 
[x] flesh out skeleton logic
[x] implement orf finder (without orfanage)
    [] create test dataset with CHESS dataset
    [] translate in 6 frames and store all orfs above threshold
    [] multithread ? 
[] integrate psauron for scoring orfs
[] implement debug in verbose mode in stderr

# Notes
- seems like pip and conda both work, package in pip, development in conda?
    - will expect the user to install orfanage, mmseqs2 themself

- orf finding steps:
    - extract long open reading frames
    - identify orfs with homology to known proteins (optional)
    - predict likely coding regions based on homology

- replicating output
    - what is the difference in output for prokaryote vs eukaryote
    

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
    - we are assuming orfs are non-overlapping, right? the current logic of the finder is maximized for efficiency using this assumption
    - what does the -c option mean? complete orfs are default?

# meeting notes w brian

- look at gc content for different organisms 
- also have a more complex start codon prediction (not a huge change overall but some are downstream)
- uses diamondblast to do the homology search
- usually go with UniRef90 as the database matching
- openBSD license?

make sure
- broad applications think of ways this could be used
- can target lots of genomes with diff characteristics
    - Plasmodium
- need to find where all the false positives are coming from
    - randomized seqeunces
- benchmark against other tools that operate in this area?
    - genemark's version of transdecoder? orfipy, getorf, orfM
- run something like glimmer and prodigal, etc.(wait they don't work they target bacteria and they need training)
- how small a peptide can you go? -> if you can demonstrate shorter peptide prediction accurately, that'd be a huge win
- the homology thing is kind of a cheat -> it's a catch-all, can fine-tune all the way, but it's not like de novo
    - want to be able to predict based on the sequence itself
    - some non-housekeeping genes, esp short peptides (signaling, defense) are not found on homology searches
- can try to properly predict things that have proteomic evidence -> data is frankly too noisy
    - try to address that in the paper, even if it doesn't work


# Links
https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=tgencodes#SG33
https://github.com/TransDecoder/TransDecoder/wiki