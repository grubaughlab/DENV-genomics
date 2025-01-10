from Bio import SeqIO
import csv
import os
import argparse
import sys

def main(sysargs = sys.argv[1:]):

    parser = argparse.ArgumentParser(add_help=False)

    parser.add_argument("--genbank-path", dest="genbank_path", help="path to general genbank submission (will end with genbank_submission)")
    parser.add_argument("--run-path", dest="run_path", help="path to consensus sequences for submission")
    parser.add_argument("-h","--help",action="store_true",dest="help")

    if len(sysargs)<1: 
        parser.print_help()
        sys.exit(0)
    else:
        args = parser.parse_args(sysargs)
        if args.help:
            parser.print_help()
            sys.exit(0)
            
    genbank_path = args.genbank_path
    run_path = args.run_path

    folder = os.path.join(genbank_path,run_path)
    unedited_aln = os.path.join(folder,"alignment.fasta")

    os.mkdir(os.path.join(folder,"qc_alignments"))

    for sero in ["DENV1", "DENV2", "DENV3", "DENV4"]:
        new_aln = open(os.path.join(folder, "qc_alignments", f"{sero}.fasta"),'w')
        for seq in SeqIO.parse(unedited_aln, "fasta"):
            if sero in seq.id:
                SeqIO.write(seq, new_aln, "fasta")
                
        new_aln.close()
    
if __name__ == '__main__':
    main()