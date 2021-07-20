'''Script to read in FASTA sequences of alignments, pull data from RCSB and SIFTS,
and then find all chain matches.'''

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

import csv, os, glob, itertools, sys, getopt
import numpy as np
from datetime import datetime
from os import path
import subprocess
import xml.etree.ElementTree as ET
import requests, json




## collect args
argString = "%s -f fasta_input_file -c chain_file"%sys.argv[0]
try:
    optlist, args = getopt.getopt(sys.argv[1:],'hf:c:')
except getopt.GetoptError:
    print('alignment2dssp -f <fasta_input_file> -c <chain_file>')
    raise Exception(argString)
    sys.exit(2)

## handle args
path = ''
chain_file = "/var/www/aln2pdb/aln2pdb/pdb_chain_ensembl.tsv"
for o, a in optlist:
    if o == '-h':
        print("="*60)
        print('usage: alignment2dssp -d <directory_of_fasta_directories> -c <chain_file>')
        print("="*60)
        sys.exit()
    if o == '-f':
        path = a
    if o == '-c':
        chain_file = a

        
def check_args():
    '''Checks the passed arguments'''
    if path == '':
        print("="*60)
        print('usage: alignment2dssp -f <fasta_input_file> -c <chain_file>')
        print("="*60)
        sys.exit()
    if chain_file == '':
        print("="*60)
        print('usage: alignment2dssp -f <fasta_input_file> -c <chain_file>')
        print("="*60)
        sys.exit()

def create_directories():
    '''Creates needed directories if they don't already exist'''
    if os.path.isdir("/var/www/aln2pdb/html/python/output/") == False:
        print("creating output directory")
        subprocess.call(["mkdir", "/var/www/aln2pdb/html/python/output"])
    if os.path.isdir("/var/www/aln2pdb/html/python/xml_files/") == False:
        print("creating SIFTS XML file directory")
        subprocess.call(["mkdir", "/var/www/aln2pdb/html/python/xml_files"])
      
        
class record_obj:
    '''Constructor for record object that will be used throughout the
    script.  The list, "input_for_aligner", is full of of these record objects.
    '''
    def __init__(self, ensembl,refseq,pdbid,chain,uniprot,pdbseq,dssp,entity,adjustment_scalar):
        self.ensembl = ensembl
        self.refseq = refseq
        self.pdbid = pdbid
        self.chain = chain
        self.uniprot = uniprot
        self.pdbseq = pdbseq
        self.dssp = dssp
        self.entity = entity
        self.adjustment_scalar = adjustment_scalar
        
    def add_bs(self,start_lig,end_lig):
        self.start_lig = start_lig
        self.end_lig = end_lig
        
    def add_ir(self,start_ir,end_ir):
        self.start_ir = start_ir
        self.end_ir = end_ir
        
    def add_mr(self,mod_res):
        self.mod_res = mod_res
    
    def __str__(self):
        print(self.ensembl, self.refseq, self.pdbid, self.chain, self.uniprot, self.pdbseq, self.dssp, self.entity, self.adjustment_scalar)
        
            
def table_builder(record): 
    '''builds output table based on multiple sequence alignments'''
    adjustment_scalar = record.adjustment_scalar
    d = []
    r = []
    local_pos = 0
    ref_pos = 0
    
    ensembl = record.ensembl
    pdb = record.pdbid
    chain = record.chain
    uniprot = record.uniprot
    ss = record.dssp
    start_lig = record.start_lig
    end_lig = record.end_lig
    start_ir = record.start_ir
    end_ir = record.end_ir
    mod_res = record.mod_res
    id_for_reference = pdb+"_"+chain
    id_for_dssp = id_for_reference+"_"+"pdb"
    reference = record.refseq
    ss_seq = record.pdbseq

    for index, residue in enumerate(reference):
        pdb_index = ref_pos-adjustment_scalar
        ligand = "-"
        interface_res = "-"
        modified_res = "-"
        if residue == "-":
            d.append(
                {
                    "ensembl_id": ensembl,
                    "uniprot_id": "-",
                    "pdb_id": "-",
                    "chain": "-",
                    "pdb_position": "-",
                    "pdb_dssp": "-",
                    "pdb_residue": residue,
                    "binding_site": ligand,
                    "interface_res": interface_res,
                    "modified_res": modified_res,
                    "fasta_position": index+1
                }

            )
        elif residue != "-":
            try:
                ligand = start_lig[local_pos]
            except:
                pass
            try:
                ligand = end_lig[local_pos]
            except:
                pass
            try:
                interface_res = start_ir[local_pos]
            except:
                pass
            try:
                interface_res = end_ir[local_pos]
            except:
                pass
            try:
                modified_res = mod_res[local_pos]
            except:
                pass
            try:
                dssp = ss[local_pos]
            except:
                dssp = "-"
            try:
                pdb_res = ss_seq[pdb_index]
            except:
                pdb_res = "-"  
                
            
            if ref_pos >= adjustment_scalar:
                d.append(
                    {
                        "ensembl_id": ensembl,
                        "uniprot_id": uniprot,
                        "pdb_id": pdb,
                        "chain": chain,
                        "ref_original_position": local_pos+1,
                        "pdb_residue": pdb_res,
                        "pdb_dssp": dssp,
                        "ref_residue": residue,
                        "binding_site": ligand,
                        "interface_res": interface_res,
                        "modified_res": modified_res,
                        "fasta_position": index+1
                    }
                )   
            
            else:
                d.append(
                    {
                        "ensembl_id": ensembl,
                        "uniprot_id": uniprot,
                        "pdb_id": pdb,
                        "chain": chain,
                        "ref_original_position": local_pos+1,
                        "pdb_residue": "-",
                        "pdb_dssp": dssp,
                        "ref_residue": residue,
                        "binding_site": ligand,
                        "interface_res": interface_res,
                        "modified_res": modified_res,
                        "fasta_position": index+1
                    }
                )
            
            ref_pos = ref_pos+1
            local_pos = local_pos+1

    df = pd.DataFrame(d)
    return df


def get_entityid(pdbid,chain):
    '''Gets the corresponding RCSB entity ID based on PDB and chain IDs'''
    r = requests.get('https://data.rcsb.org/rest/v1/core/polymer_entity_instance/'+pdbid+'/'+chain)
    results = r.json()
    try:
        entityid = results["rcsb_polymer_entity_instance_container_identifiers"]["entity_id"]
        adjustment_scalar = int(results["rcsb_polymer_entity_instance_container_identifiers"]["auth_to_entity_poly_seq_mapping"][0])-1
    except: # No entity ID found. Check to see if the aysm_id is the same as auth_asym_id
        entityid = None
        adjustment_scalar = None

    return entityid, adjustment_scalar


def get_dssp(pdbid,chain):
    '''Gets secondary structure information and PDB sequence using RCSB's new REST API '''
    ss = {}
    pdb_seq = ''
    entityid,adjustment_scalar = get_entityid(pdbid,chain)
    
    if entityid == None:
        return "","",entityid,adjustment_scalar
    
    else:
        r = requests.get('https://data.rcsb.org/rest/v1/core/polymer_entity/'+pdbid+'/'+entityid)
        results = r.json()
        pdb_seq = results["entity_poly"]["pdbx_seq_one_letter_code_can"]

        valid = ["loop","strand","sheet","helix","turn"]
        s = requests.get('https://data.rcsb.org/rest/v1/core/polymer_entity_instance/'+pdbid+'/'+chain)
        dssp_results = s.json()

        for i in dssp_results["rcsb_polymer_instance_feature"]:
            try:
                ss_name = i['name']
                if ss_name in valid:
                    for feature in i['feature_positions']:
                        ss_start = feature['beg_seq_id']
                        ss_end = feature['end_seq_id']
                        ss_start = ss_start+adjustment_scalar
                        ss_end = ss_end+adjustment_scalar
                        for x in range(ss_start-1, ss_end):
                            ss[x] = ss_name
            except:
                pass

        return pdb_seq, ss, entityid, adjustment_scalar


def fasta_parser(filename):
    '''Parses the multiple alignment fasta file and returns a unique list of ensembl ids, pdb ids,
    uniprot ids, chains, and sequences.  Requires a PDB chain file to match ensembl IDs.
    input_for_aligner is full of record objects and are:
    -Ensembl ID
    -reference sequence (from fasta file)
    -matched PDB ID from chain file
    -matched chain ID from chain file
    -matched UniProt ID from chain file
    -RCSB sequence from matched PDB ID 
    -RCSB DSSP sequence that corresponds to RCSB sequence from matched PDB ID
    '''
    print("opening",filename)
    records = []
    chains = []
    input_for_aligner = []
    ss_seq = ""
    dssp = ""
    entityid = ""

    for index, record in enumerate(SeqIO.parse(filename, "fasta")):
        seq = str(record.seq)
        ensembl_id = record.id
        with open(chain_file, 'r') as f:
            for line in f.readlines():
                if ensembl_id in line:
                    out = line.split("\t")
                    chains.append([ensembl_id,seq,out[0],out[1],out[2]])

    b_set = set(map(tuple,chains))  #need to convert the inner lists to tuples so they are hashable
    unique_chains = list(b_set)
    print("number of unique chains:",len(unique_chains))
    if len(unique_chains) == 0:
        return []
    
    else:
        for chain in unique_chains:
            pdb_id = chain[2]
            chain_X = chain[3].upper()
            ss_seq = ''
            dssp = ''
            ss_seq_id = pdb_id.upper()+":"+chain_X+":"+"sequence"
            ss_str_id = pdb_id.upper()+":"+chain_X+":"+"secstr"
            ss_seq, dssp, entityid, adjustment_scalar = get_dssp(pdb_id,chain_X)
            
            if len(ss_seq) == 0: 
                print("No PDB sequence found for %s / %s. Check that aysm_id is the same as auth_asym_id" % (pdb_id,chain_X))
                continue
            if len(dssp) == 0:
                print("No DSSP sequence found for %s / %s. Check that aysm_id is the same as auth_asym_id" % (pdb_id,chain_X))
                dssp = "-"
            record_object = record_obj(chain[0],chain[1],chain[2],chain[3],chain[4],ss_seq,dssp,entityid,adjustment_scalar)
            input_for_aligner.append(record_object)

        return input_for_aligner


def sifts_parse(record):
    '''Downloads and parses multiple SIFTS XML records for a given set of parsed fasta records.
    Returns structural info from the SIFTS record.'''
    d = []
    protein = record.pdbid
    target_chain = record.chain
    parent_dir = protein[1:3]
    ftp_path = "ftp.ebi.ac.uk/pub/databases/msd/sifts/split_xml/"
    gz_file = protein+".xml.gz"
    file_path = ftp_path+parent_dir+"/"+gz_file
    xml_file = protein+".xml"
    
    if os.path.isfile("/var/www/aln2pdb/html/python/xml_files/"+xml_file) == False:
        try:
            print("downloading",protein,target_chain,"from SIFTS")
            # download SIFTS xml file
            subprocess.call(["wget", file_path])
            subprocess.call(["gzip","-d", gz_file])
            subprocess.call(["mv", xml_file,"/var/www/aln2pdb/html/python/xml_files/"])
        except:
            pass
    else:
        print("loading",protein,target_chain,"from SIFTS")

    # parse XML
    tree = ET.parse("/var/www/aln2pdb/html/python/xml_files/"+xml_file)
    root = tree.getroot()
    
    # finding the correct entity for the chain
    found_chain = "unknown"
    for entity in root.iter('{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}entity'):
        entityid = entity.attrib['entityId'] 
        for residue in root.findall("./{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}entity[@entityId='"+entityid+"']/{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}segment/{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}listResidue/{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}residue"):
            for details in residue.findall('./{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}residueDetail/[@property="nameSecondaryStructure"]'):
                structure = details.text
            found_chain_element = residue.find('./{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}crossRefDb')
            found_chain = found_chain_element.attrib['dbChainId']

            if found_chain == target_chain:
                found_entityid = entityid #SIFTS uses letters for entity IDs, and might not correspond to chain ID
                break
    
    for residue in root.findall("./{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}entity[@entityId='"+found_entityid+"']/{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}segment/{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}listResidue/{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}residue"):
        for crossDb in residue.findall("./{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}crossRefDb[@dbSource='PDB']"):
            pass
        for details in residue.findall('./{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}residueDetail/[@property="nameSecondaryStructure"]'):
            structure = details.text

        if crossDb.attrib['dbResNum'] == 'null':
            continue
        else:
            try:
                d.append(
                    {
                        "authPos": int(crossDb.attrib['dbResNum']),
                        "residue": crossDb.attrib['dbResName'],
                        "sifts_structure": structure
                    }

                )
            except Exception as e:
                try:
                    d.append(
                        {
                            "authPos": int(crossDb.attrib['dbResNum']),
                            "residue": crossDb.attrib['dbResName'],
                            "sifts_structure": "-"
                        }

                    )
                except:
                    pass

    df = pd.DataFrame(d)
    return df


def get_best_ss(results):
    '''uses a mode function to find the most frequent result for sifts, pdb, and original position
    from the results. Accounts for two-way ties'''
    df = results
    pdb = df[["pdb_dssp"]]
    sifts = df[["sifts_structure"]]
    position = df[["fasta_position"]]

    pdb.replace(r'^\s+$', np.nan, regex=True,inplace=True)
    pdb.replace('-', np.nan, regex=True,inplace=True)
    sifts.replace('-', np.nan, regex=True,inplace=True)
    position.replace('-', np.nan, regex=True,inplace=True)
    
    best_sifts = sifts.mode(axis=1)
    best_pdb = pdb.mode(axis=1)
    best_position = position.mode(axis=1)

    # accounting for two-way ties
    if best_sifts.shape[1] == 2:
        best_sifts.columns=["top_sifts_structure","alt_top_sifts_structure"]
    else:
        best_sifts=best_sifts.rename(columns = {0:'top_sifts_structure'})
    if best_pdb.shape[1] == 2:
        best_pdb.columns=["top_pdb_dssp","alt_top_pdb_dssp"]
    else:
        best_pdb=best_pdb.rename(columns = {0:'top_pdb_dssp'})
    if best_position.shape[1] == 2:
        best_position.columns=["original_position","alt_position"]
    else:
        best_position=best_position.rename(columns = {0:'original_position'})

    best_pdb.replace(np.nan, "-",inplace=True)
    best_sifts.replace(np.nan, "-",inplace=True)
    best_position.replace(np.nan, "-",inplace=True)

    summarized_results = pd.concat([best_position, best_pdb, best_sifts], axis=1)
    return summarized_results


def binding_sites(input_for_aligner):
    '''Queries the PBDe API for binding sites'''
    for record in input_for_aligner:
        start_lig = {}
        end_lig = {}
        pdbid = record.pdbid
        entityid = record.entity
        adjustment_scalar = record.adjustment_scalar

        r = requests.get('https://www.ebi.ac.uk/pdbe/graph-api/pdbe_pages/binding_sites/'+pdbid+'/'+str(entityid))

        response = r.json()
        if bool(response) == False:
            record.add_bs(start_lig,end_lig) # empty dictionaries for start and end sites
            continue

        pdb_seq = record.refseq

        for ligand in response[pdbid]["data"]:
            lig_name = ligand['name']
            for i in ligand['residues']:
                startI, startC = i['startIndex'],i['startCode']
                endI, endC = i['endIndex'],i['endCode']
                startI = startI+adjustment_scalar
                endI = endI+adjustment_scalar
                try:
                    start_lig[startI-1] = lig_name
                    end_lig[endI-1] = lig_name
                except Exception as e:
                    pass

        record.add_bs(start_lig,end_lig)


def interface_residues(input_for_aligner):
    '''Queries the PBDe API for interface residues'''
    for record in input_for_aligner:
        start_ir = {}
        end_ir = {}
        pdbid = record.pdbid
        entityid = record.entity
        adjustment_scalar = record.adjustment_scalar

        r = requests.get('https://www.ebi.ac.uk/pdbe/graph-api/pdbe_pages/interfaces/'+pdbid+'/'+str(entityid))

        response = r.json()
        if bool(response) == False:
            record.add_ir(start_ir,end_ir) # empty dictionaries for start and end sites
            continue

        pdb_seq = record.refseq

        for ir in response[pdbid]["data"]:
            ir_name = ir['name']
            for i in ir['residues']:
                startI, startC = i['startIndex'],i['startCode']
                endI, endC = i['endIndex'],i['endCode']
                startI = startI+adjustment_scalar
                endI = endI+adjustment_scalar
                try:
                    start_ir[startI-1] = ir_name
                    end_ir[endI-1] = ir_name
                except Exception as e:
                    pass

        record.add_ir(start_ir,end_ir)


def modified_residues(input_for_aligner):
    '''Queries the PBDe API for modified residues'''
    for record in input_for_aligner:
        mod_res = {}
        pdbid = record.pdbid
        chain = record.chain
        adjustment_scalar = record.adjustment_scalar

        r = requests.get('https://www.ebi.ac.uk/pdbe/graph-api/pdb/modified_AA_or_NA/'+pdbid)
        response = r.json()
        if bool(response) == False:
            record.add_mr(mod_res) # empty dictionaries for start and end sites
            continue

        pdb_seq = record.refseq
        
        for i in (response[pdbid]):
            if i["chain_id"] == chain:
                chem_name = i["chem_comp_name"]
                local_pos = i["residue_number"]
                local_pos = local_pos+adjustment_scalar
                try:
                    pdb_seq[local_pos-1]
                    mod_res[local_pos-1] = chem_name
                except Exception as e:
                    pass

        record.add_mr(mod_res)


if __name__ == "__main__": 
    startTime = datetime.now()
    success = 0
    fail = 0
    #check_args()
    #create_directories()
    #fasta = "/var/www/aln2pdb/aln2pdb/fastas/ENSGT00390000005995_1.fasta.fas"
    fasta = path
    print(fasta)
    results = pd.DataFrame({'A' : []}) # create empty dataframe
    id_mapper = {}
    loopTime = datetime.now()
    filename = fasta
    input_for_aligner = fasta_parser(filename)
    if len(input_for_aligner) == 0:
        print("No pdb chain hits for",filename)
        print("Execution time:",(datetime.now() - startTime))
        print("="*60)
        fail = fail+1
    print("Getting binding sites...")
    binding_sites(input_for_aligner)
    print("Getting interface residues...")
    interface_residues(input_for_aligner)
    print("Getting modified residues...")
    modified_residues(input_for_aligner)
    for record in input_for_aligner:
        found_entityid = ""
        if results.empty == True:
            results = table_builder(record)
            sifts = sifts_parse(record)
            results = pd.merge(results,sifts,how='left',left_on="ref_original_position",right_on="authPos")
            results = results.replace(to_replace = np.nan, value ="-") 
        else:
            tmp = table_builder(record)
            sifts = sifts_parse(record)
            tmp = pd.merge(tmp,sifts,how='left',left_on="ref_original_position",right_on="authPos")
            tmp = tmp.replace(to_replace = np.nan, value ="-") 
            results = pd.concat([results,tmp],axis=1)
        #hits.append(record.ensembl+":"+record.pdbid+":"+record.chain)

    results.replace("-",np.nan, regex=True,inplace=True)
    results.to_csv("/var/www/aln2pdb/html/python/output/"+filename.split("/")[-1].split(".fas")[0]+"_dssp.txt",sep="\t",index=False)
    summarized_results=get_best_ss(results)
    summarized_results.to_csv("/var/www/aln2pdb/html/python/output/"+filename.split("/")[-1].split(".fas")[0]+"_summarized.txt",sep="\t",index=False)

    print("-"*60)
    #print("Number of hits to PBD:",len(hits))
    #print(", ".join(hits))
    #print("Execution time:",(datetime.now() - loopTime))
    #print("="*60)
    #success = success+1

            
    #print("successful runs:",success)
    #print("null runs:",fail)
    print("Total execution time:",(datetime.now() - startTime))
    
    
