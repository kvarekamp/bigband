import dataclasses
from pathlib import Path
import logging
import gzip
import types
from typing import List, Text
from datetime import datetime
from dateutil.parser import parse
from dataclasses import dataclass
from typing import TextIO
import json

log = logging.getLogger('genbank_reader')
log.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(logging.Formatter('%(asctime)s|%(name)s|%(levelname)s| %(message)s'))
log.addHandler(ch)

@dataclass
class GenbankHeader():
    file_name: str
    datetime: datetime
    version_major: int
    version_minor: int
    n_loci: int
    n_bases: int
    n_sequences: int

# With the exception of the lists of new, changed, and
# deleted accession numbers, each of the files of a GenBank release begins
# with the same header, except for the first line, which contains the file
# name, and the sixth line, which contains the title of the file. The first
# line of the file contains the file name in character positions 1 to 9 and
# the full database name (Genetic Sequence Data Bank, aka 'GenBank') starting
# in column 22. The brief names of the files in this release are listed in
# Section 2.2.

#   The second line contains the date of the current release in the form
# `day month year', beginning in position 27. The fourth line contains
# the current GenBank release number. The release number appears in
# positions 48 to 52 and consists of three numbers separated by a decimal
# point. The number to the left of the decimal is the major release
# number. The digit to the right of the decimal indicates the version of
# the major release; it is zero for the first version. The sixth line
# contains a title for the file. The eighth line lists the number of
# entries (loci), number of bases (or base pairs), and number of reports
# of sequences (equal to number of entries in this case). These numbers are
# right-justified at fixed positions. The number of entries appears in
# positions 1 to 8, the number of bases in positions 16 to 26, and the
# number of reports in positions 40 to 47. The third, fifth, seventh, and
# ninth lines are blank.
def read_genbank_header(r: TextIO) -> GenbankHeader:    
    lines = []
    for i in range(9):
        lines.append(next(r))

    version_major, _, version_minor = lines[3][47:].strip().partition(".")
    str_loci = lines[7][:8].strip()
    str_bases = lines[7][15:26].strip()
    str_sequences = lines[7][39:47].strip()
    return GenbankHeader(
        file_name=lines[0][:9].strip(),
        datetime=parse(lines[1]),
        version_major=int(version_major),
        version_minor=int(version_minor),
        n_loci=int(str_loci),
        n_bases=int(str_bases),
        n_sequences=int(str_sequences),
    )


def read_genbank_seqs(genbank:Path):
    one_dumped = False

    for p in genbank.glob("*.seq.gz"):
        seqs = read_genbank_seq_file(p)

        if not one_dumped:
            with open('dump.json', 'w') as fout:
                json.dump(seqs, fout, indent=2)
                one_dumped = True


def read_genbank_seq_file(p:Path, cleanup_origin=True):
    """Reads a genbank file according to the file format desribed in 
    https://www.ncbi.nlm.nih.gov/genbank/release/current/ """

    with gzip.open(p,'rt') as fin:
        h = read_genbank_header(fin)
        loci = []
        locus = {}
        last_h = None
        last_h1 = None
        refcount = 0
        for line in fin:
            #first work out the line type and colwidth
            if line.isspace():
                continue
            elif line[:12].isspace():
                ltype = "TXT"
            elif line[:1].isspace():
                ltype = "H2"
            elif line[:2] == "//":
                ltype = "END"
            else:
                ltype = "H1"
                if line.startswith("FEATURES"):
                    colwidth = 21
                elif line.startswith("ORIGIN"):
                    colwidth = 9
                else:
                    colwidth = 12

            #get the header and do some cleanup on it
            header = line[:colwidth].strip()
            if header.isdigit():
                header = line[:colwidth].replace(' ', '0')
            elif header == "REFERENCE":
                refcount += 1
                header = f"REFERENCE.{refcount}"

            content = line[colwidth:].strip()
            if ltype == "END":
                loci.append(locus)
                locus = {}
                last_h = None
                last_h1 = None
                refcount = 0
            if ltype == "TXT":
                locus[last_h] += "\n"+content
            elif ltype == "H2":
                last_h = f"{last_h1}.{header}"
                locus[last_h] = content
            elif ltype == "H1":
                last_h1 = header 
                last_h = header 
                locus[last_h] = content

        if cleanup_origin:
            for locus in loci:
                origins = sorted(key for key in locus.keys() if key.startswith("ORIGIN."))
                # log.debug(f"cleaning up {len(origins)} origins")
                lines = [locus[origin] for origin in origins]
                locus["ORIGIN"] = "".join(lines).replace(" ", "")
                for origin in origins:
                    del locus[origin]

        log.debug(f"read {len(loci)} loci from {p}")
        return loci
    


if __name__ == "__main__":
    read_genbank_seqs(Path("genbank"))
