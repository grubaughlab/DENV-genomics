import os
import sys
import csv
from collections import defaultdict
import datetime as dt
from unidecode import unidecode
from Bio import SeqIO
import tqdm as tqdm
import argparse

def main(sysargs = sys.argv[1:]):

    parser = argparse.ArgumentParser(add_help=False)

    parser.add_argument("--master-csv", dest="master_csv", help="GLab master csv")
    parser.add_argument("--file-path", dest="file_path", help="Overall path to where genomic data is kept, relative to script")
    parser.add_argument("--genbank-path", dest="genbank_path", help="path to where genbank files will be made, relative to current directory")
    parser.add_argument("--date-submission", dest="date_submission", help="for file name within lab folder")
    parser.add_argument("-h","--help",action="store_true",dest="help")

    if len(sysargs)<1: 
        parser.print_help()
        sys.exit(0)
    else:
        args = parser.parse_args(sysargs)
        if args.help:
            parser.print_help()
            sys.exit(0)

    master_csv = args.master_csv
    file_path = args.file_path
    genbank_path = args.genbank_path
    date_submission = args.date_submission


    run_to_seqs, new_names, seqs_to_lab, lab_to_seqs, all_labs, seqs_to_fasta_name, seq_metadata = parse_master_csv(master_csv)
    full_genbank_paths = prep_folders(all_labs, genbank_path, date_submission)
    
    pull_fasta_files(all_labs, run_to_seqs, full_genbank_paths, file_path, seqs_to_lab, seqs_to_fasta_name, new_names)

    make_sra_metadata(lab_to_seqs, full_genbank_paths, seq_metadata, new_names)
    make_biosample_metadata(seq_metadata, new_names, lab_to_seqs, full_genbank_paths)
    make_genbank_metadata(seq_metadata, new_names, lab_to_seqs, full_genbank_paths)


    
def parse_master_csv(master_csv):
    
    #parse the metadata to get everything that's needed
    run_to_seqs = defaultdict(list)
    new_names = {}
    seqs_to_lab = {}
    lab_to_seqs = defaultdict(list)
    seqs_to_fasta_name = {}
    seq_metadata = defaultdict(dict)

    with open(master_csv) as f:
        data = csv.DictReader(f)
        for l in data:
            if l['Final_Genome'] == "1" and l['GenBank_Accession'] == "":
                
                name = l['Yale-ID']

                if l['NGS_Serotype'] == "":
                    print(f"Genome marked as best but no serotype from sequencing: {name}\n")
                    continue
                
                if "+" in l['NGS_Run_ID']:
                    run = l['NGS_Run_ID'].split("+")[1].lstrip(" ")
                else:
                    run = l['NGS_Run_ID']
                            
                run_to_seqs[run].append(name)
                
                
                if l['Country'] == "USA":
                    location = f'USA:{l["Division_(state)"]}'
                else:
                    location = l["Country"].replace("; ", "-").replace(", ", "-").replace(" ","_")

                if ";" in location or "," in location:
                    print(f"Multiple locations detected for {name} in {l['NGS_Serotype']}: please check and edit CSV if necessary")

                new_name = f'{l["NGS_Serotype"]}/{location}/{l["Collection_Date"]}/{name}'
                if len(new_name) > 50:
                    print("fasta names are too long for genbank")
                new_names[name] = new_name
                
                lab = unidecode(l["Lab_Source"].replace(" ","_"))
                seqs_to_lab[name] = lab
                lab_to_seqs[lab].append(name)
                
                if "_P" in run:
                    new_id = "-".join(name.split("-")[0:2])
                    seqs_to_fasta_name[name] = f'{new_id}P.{l["NGS_Serotype"]}.20.cons.fa'
                    bam_name = f'{new_id}P.bam'
                else:
                    seqs_to_fasta_name[name] = f'{name}.{l["NGS_Serotype"]}.20.cons.fa'
                    bam_name = f'{name}.bam'

                
                seq_metadata[name]["bam_name"] = bam_name
                seq_metadata[name]["date"] = l['Collection_Date']
                seq_metadata[name]["location"] = location.replace("_"," ")
                seq_metadata[name]["serotype"] = l['NGS_Serotype']
                seq_metadata[name]["genotype"] = l['NGS_Genotype']
                
                if l['Location_diagnosed'].replace(" ","_") != location:
                    if "-" in location or ";" in location:
                        location_caught = location.replace("-", " or ").replace(";", " or ")
                    else:
                        location_caught = location.replace('_', ' ')
                    diagnosed = f"traveller from {location_caught} diagnosed in {l['Location_diagnosed'].replace('_', ' ')}"
                else:
                    diagnosed = ""
                
                seq_metadata[name]["diagnosed"] = diagnosed
                seq_metadata[name]["original_id"] = l['Original_ID']
                
            if "/" in l['Collection_Date']:
                print("CSV needs reformatting for dates")
                break
                
        all_labs = set(list(seqs_to_lab.values()))
        print(f'{len(seq_metadata)} prepped for submission')

        return run_to_seqs, new_names, seqs_to_lab, lab_to_seqs, all_labs, seqs_to_fasta_name, seq_metadata

def prep_folders(all_labs, genbank_path, date_submission):
    ## prep folders
    full_genbank_paths = {}
    for i in all_labs:
        if i not in os.listdir(genbank_path):
            os.mkdir(os.path.join(genbank_path,i))

        if date_submission not in os.listdir(os.path.join(genbank_path,i)):
            os.mkdir(os.path.join(genbank_path,i,date_submission))
            
        full_genbank_paths[i] = os.path.join(genbank_path,i,date_submission)

    return full_genbank_paths

def pull_fasta_files(all_labs, run_to_seqs, full_genbank_paths, file_path, seqs_to_lab, seqs_to_fasta_name, new_names):
    lab_fasta = {}
    for lab in all_labs:
        lab_fasta[lab] = open(os.path.join(full_genbank_paths[lab], "alignment.fasta"), 'w')

    print("pulling samples into fasta files")
    for run, sequences in tqdm.tqdm(run_to_seqs.items()):
        if "_P" in run:
            full_file_path = os.path.join(file_path, "Pan-serotype", run, "CONSENSUS")
        else:
            full_file_path = os.path.join(file_path, run, "CONSENSUS")
            
        for seq in sequences:
            fasta_file = lab_fasta[seqs_to_lab[seq]]
            seq_file = os.path.join(full_file_path, seqs_to_fasta_name[seq])
            for seq_obj in SeqIO.parse(seq_file, "fasta"):
                seq_obj.name = new_names[seq]
                seq_obj.id = seq_obj.name
                seq_obj.description = seq_obj.name
                SeqIO.write(seq_obj, lab_fasta[seqs_to_lab[seq]], "fasta")
                  
    for lab, fasta in lab_fasta.items():
        fasta.close()

    return

def make_sra_metadata(lab_to_seqs, full_genbank_paths, seq_metadata, new_names):

    print("making sra metadata")

    ref_genome = {"DENV1":"NC_001477.1","DENV2":"NC_001474.2","DENV3":"NC_001475.2","DENV4":"NC_002640.1"}

    headers = ["sample_name", "library_ID", "title", "library_strategy", "library_source", "library_selection", "library_layout",
            "platform", "instrument_model", "design_description", "filetype", "filename", "filename2", "filename3", "filename4",
            "assembly", "fasta_file"]

    for lab, sequences in lab_to_seqs.items():    
        with open(os.path.join(full_genbank_paths[lab], "sra_metadata.tsv"), 'w') as fw:
            writer = csv.DictWriter(fw, fieldnames=headers, delimiter="\t")
            writer.writeheader()
            
            for seq in sequences:
                write_dict = {}
                
                write_dict["sample_name"] = new_names[seq]
                write_dict["library_ID"] = new_names[seq]
                write_dict["title"] = "Yale dengue virus genomic epidemiology"
                write_dict["library_strategy"] = "AMPLICON"
                write_dict["library_source"] = "VIRAL RNA"
                write_dict["library_selection"] = "PCR"
                write_dict["library_layout"] = "paired"
                write_dict["platform"] = "ILLUMINA"
                write_dict["instrument_model"] = "Illumina NovaSeq 6000"
                write_dict["design_description"] = "Amplicon-based whole genome sequencing of dengue virus"
                write_dict["filetype"] = "bam"
                write_dict["filename"] = seq_metadata[seq]["bam_name"]
                write_dict["assembly"] = ref_genome[seq_metadata[seq]["serotype"]]
                
                for i in headers:
                    if i not in write_dict:
                        write_dict[i] = ""
                        
                writer.writerow(write_dict)

    
def make_biosample_metadata(seq_metadata, new_names, lab_to_seqs, full_genbank_paths):

    print("making biosample metadata")

    headers = ["sample_name", "sample_title", "bioproject_accession", "organism", "strain", "isolate", "collected_by",
        "collection_date", "geo_loc_name", "host", "host_disease", "isolation_source", "lat_lon", "culture_collection",
        "genotype", "host_age", "host_description", "host_disease_outcome", "host_disease_stage", "host_health_state",
        "host_sex", "host_subject_id", "host_tissue_sampled", "passage_history", "pathotype", "serotype", "serovar",
        "specimen_voucher", "subgroup", "subtype", "description"]

    for lab, sequences in lab_to_seqs.items():    
        with open(os.path.join(full_genbank_paths[lab], "biosample_metadata.tsv"), 'w') as fw:
            writer = csv.DictWriter(fw, fieldnames=headers, delimiter="\t")
            writer.writeheader()
            
            for seq in sequences:
                write_dict = {}
                
                write_dict["sample_name"] = new_names[seq]
                write_dict["organism"] = "Dengue virus"
                write_dict["isolate"] = seq
                write_dict["collected_by"] = lab
                write_dict["collection_date"] = seq_metadata[seq]["date"]
                write_dict["geo_loc_name"] = seq_metadata[seq]["location"]
                write_dict["host"] = "Homo sapiens"
                write_dict["host_disease"] = "Dengue"
                write_dict["isolation_source"] = "Serum"
                write_dict["lat_lon"] = "Unknown"
                write_dict["genotype"] = seq_metadata[seq]["genotype"]
                write_dict["serotype"] = seq_metadata[seq]["serotype"][-1]
                write_dict["description"] = seq_metadata[seq]["diagnosed"]

                for i in headers:
                    if i not in write_dict:
                        write_dict[i] = ""
                        
                writer.writerow(write_dict)


def make_genbank_metadata(seq_metadata, new_names, lab_to_seqs, full_genbank_paths):

    print("making genbank metadata")

    headers = ["Sequence_ID", "isolate", "country", "host", "collection-date", "isolation-source", "serotype", "genotype",
          "strain", "notes"]

    for lab, sequences in lab_to_seqs.items():    
        with open(os.path.join(full_genbank_paths[lab], "genbank_source_modifiers.tsv"), 'w') as fw:
            writer = csv.DictWriter(fw, fieldnames=headers, delimiter="\t")
            writer.writeheader()
            
            for seq in sequences:
                write_dict = {}
                
                write_dict["Sequence_ID"] = new_names[seq]
                write_dict["isolate"] = new_names[seq]
                write_dict["country"] = seq_metadata[seq]["location"]
                write_dict["host"] = "Homo sapiens"
                write_dict["collection-date"] = seq_metadata[seq]["date"]
                write_dict["isolation-source"] = "Serum"
                write_dict["serotype"] = f'{seq_metadata[seq]["genotype"]} genotype'
                write_dict["genotype"] = seq_metadata[seq]["serotype"][-1]
                write_dict["strain"] = seq_metadata[seq]["original_id"]
                write_dict["notes"] = seq_metadata[seq]["diagnosed"]
                
                
                writer.writerow(write_dict)
            

    
if __name__ == '__main__':
    main()