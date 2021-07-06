from pathlib import Path
from fasta_funcs import read_fastas
from tokenizers import Tokenizer
from tokenizers.models import WordPiece
from tokenizers.trainers import WordPieceTrainer


tokenizer = Tokenizer(WordPiece(unk_token="-"))
trainer = WordPieceTrainer(special_tokens=["-"])
tokenizer.train_from_iterator(read_fastas(Path("genomes"), max_lines=5e6), trainer)
tokenizer.save("tokenizer-bigband.json")
