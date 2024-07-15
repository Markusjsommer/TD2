conda activate td
pip install .

TD2.LongOrfs -t ./tests/data/MANE_test_rna.fna --m-start -O ./results/
python ./tests/compare_peps.py ./tests/data/MANE_test_rna.fna.transdecoder_dir/longest_orfs.pep ./results/transcripts.transmark_dir/longest_orfs.pep --diff --matched --internal