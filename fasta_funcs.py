from ftplib import FTP
from pathlib import Path
from collections import Counter
import gzip
# from the readme:
#   *_genomic.fna.gz file
#        FASTA format of the genomic sequence(s) in the assembly. Repetitive 
#        sequences in eukaryotes are masked to lower-case (see below).
#        The FASTA title is formatted as sequence accession.version plus 
#        description. The genomic.fna.gz file includes all top-level sequences in
#        the assembly (chromosomes, plasmids, organelles, unlocalized scaffolds,
#        unplaced scaffolds, and any alternate loci or patch scaffolds). Scaffolds
#        that are part of the chromosomes are not included because they are
#        redundant with the chromosome sequences; sequences for these placed 
#        scaffolds are provided under the assembly_structure directory.

def download_fastas(dir):
    ftp = FTP("ftp.ncbi.nlm.nih.gov", timeout=3600)
    ftp.login()
    for remotedir in ('genomes/all/GCA', 'genomes/all/GCF'):
        for subdir in ftp.nlst(remotedir):            
            for subsubdir in ftp.nlst(subdir):
                if subsubdir < "genomes/all/GCA/000/011":
                    print(f"skipping subsubdir {subsubdir}")
                    continue
                for subsubsubdir in ftp.nlst(subsubdir):
                    for assemblydir in ftp.nlst(subsubsubdir):
                        try: 
                            items = ftp.nlst(assemblydir)
                        except Exception as e:
                            print(f"could not list {assemblydir}: {e}")
                            continue
                        for item in items:
                            if item.endswith("_genomic.fna.gz") and not item.endswith("_from_genomic.fna.gz"):
                                try:
                                    print(item)
                                    fpath = dir / Path(item).name
                                    if not fpath.exists():
                                        print(f"downloading to {fpath}")
                                        with fpath.open('wb') as fout:
                                            ftp.retrbinary(f"RETR {item}", fout.write)
                                    else:
                                        print(f"I already have {fpath}, skipping")
                                except Exception as e:
                                    print(f"could not download {item}: {e}")

# example headers
# >lcl|AE014298.5_cds_ABW09315.1_1 [gene=CG17636] [locus_tag=Dmel_CG17636] [db_xref=FLYBASE:FBpp0111834] [protein=uncharacterized protein, isoform A] [protein_id=ABW09315.1] [location=complement(join(124464..125409,125495..126259,126626..126630))] [gbkey=CDS]
# >lcl|AE014298.5_cds_AFH07158.1_2 [gene=CG17636] [locus_tag=Dmel_CG17636] [db_xref=FLYBASE:FBpp0297938] [protein=uncharacterized protein, isoform C] [protein_id=AFH07158.1] [location=complement(join(124464..125409,125495..126259,126363..126409))] [gbkey=CDS]
def read_fastas(dir, max_lines=-1):
    nlines = 0
    for i, fpath in enumerate(dir.glob("*_genomic.fna.gz")):
        try:
            with gzip.open(fpath,'rt') as fin:
                for line in fin:
                    line = line.strip()
                    if line.startswith(">"):
                        continue
                    if max_lines > -1 and nlines >= max_lines:
                        return
                    nlines += 1
                    yield line
        except Exception as e:
            print(f"could not read {fpath}: {e}")



if __name__ == "__main__":
    download_fastas(Path("genomes"))
