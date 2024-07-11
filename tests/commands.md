conda activate td
pip install .

transmark.LongOrfs -t ./data/MANE_smallset_rna.fna --m-start -O ./results/
python ./tests/compare_peps.py ./data/MANE_smallset_rna.fna.transdecoder_dir/longest_orfs.pep ./results/transcripts.transmark_dir/longest_orfs.pep --diff --matched --internal