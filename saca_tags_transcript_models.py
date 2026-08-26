import pysam
import pandas as pd
import glob
import os
import gzip
import argparse
import sys
from functools import partial
from concurrent.futures import ProcessPoolExecutor
from bisect import bisect_right
from collections import defaultdict

# A futuro esto deberíamos cambiarlo desde el principio o dar tu un dic
PROB_LIM=0.95
cod_mod_largo = {"m":"m5C","a":"m6A", "19228":"2OmeC", "17802":"pseU","19227":"2OmeU", "17596":"inosine" , "69426":"2OmeA"}
nombre_a_codigo = {v: k for k, v in cod_mod_largo.items()}
cod_mod_1 = {"m":"m","a":"a", "19228":"C", "17802":"P","19227":"U", "17596":"I" , "69426":"A"}
ASCII_MAP = {p: chr(33 + round((p * 93) / 256)) for p in range(257)}

def prob_ASCII(prob:float) -> str:
    '''
    We convert the prob of the mods (from 0 to 255) to a character 
    in the range of 33 to 126 (printable ASCII characters)
    '''
    if prob >= 0 and prob <= 256:
        ascii_code = chr(33 + round((prob * 93) / 256))
    else:
        print(prob)
        raise ValueError
    return ascii_code

def cargaBed(path: str) -> set:
    '''
    We load the bed file with filtered sites and return a set of tuples 
    (chrom, start_position, strand)
    '''
    conf_sites = set()
    with open(path, "r") as f:
        for line in f:
            # Si la linea es commentario o está vacia
            if line.startswith("#") or not line.strip():
                continue
            cols = line.strip().split("\t")
            # Añadimos una tupla (chrom,start_position,strand)
            conf_sites.add((cols[0],int(cols[1]),cols[5]))
    return conf_sites

def cargaBed(path: str) -> dict:
    '''
    We load the bed file with filtered sites and return a dictionary 
    keyed by (chrom, strand) containing a set of positions.
    '''
    conf_sites = defaultdict(set)
    with open(path, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.strip().split("\t")
            conf_sites[(cols[0], cols[5])].add(int(cols[1]))
    return conf_sites

def cargaRead_Transcrito(path:str) -> dict:
    '''
    We extract the read_id to transcript_id association from the relations 
    file. It recieves a path to a tsv file and returns a dictionary 
    {read_id: transcript_id}
    '''
    read_tx = {}
    open_fn = gzip.open if path.endswith("gz") else open
    with open_fn(path,"rt") as f:
        for line in f:
            if line.startswith("#") or line.startswith("read_id"):
                continue
            cols = line.strip().split("\t")
            if len(cols) == 2 and cols[1] != "*":
                read_tx[cols[0]] = cols[1]
    return read_tx

def cargaGTF(gtf_path: str) -> dict:
    '''
    We extract from the gtf file the transcript structure and return a dict
    with the structure:
    {transcript_id: 
        {"chrom": chrom, 
        "strand": strand, 
        "exons": [(ini, fin), ...]}, 
        "exons_offset": [(ini, fin, offset), ...], 
        "total_length": total_length}} 
    where exons are sorted by strand and exons_offset contains the 
    length of exons before the current one, and total_length is the total 
    length of the transcript.
    '''
    tx_dict = {}
    with open(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            # Col 0: chr
            # Col 2: tipo (transcript/exon)
            # Col 3: Ini (1-based)
            # Col 4: fin (1-based e intervalo cerrado --> pos fin se incluye)
            # Col 6: Strand
            # Col 8: gene_id "ENSGXXX"; transcript_id "ENSTXXX; exon_number "X"; exon_id "ENSEXXX"
            cols = line.strip().split("\t")
            chrom = cols[0]
            # Pasamos a 0 based con pos fin no incluida
            ini = int(cols[3]) - 1
            fin = int(cols[4])
            strand = cols[6]
            # Para extraer el transcript id de aqui
            attr_string = cols[8]
            tx_id = None
            for attr in attr_string.split(";"):
                attr = attr.strip()
                if attr.startswith("transcript_id"):
                    tx_id = attr.split('"')[1]
                    break
            if not tx_id:
                continue
            if tx_id  not in tx_dict.keys():
                tx_dict[tx_id] = {"chrom":chrom,"strand":strand,"exons":[]}
            tx_dict[tx_id]["exons"].append((ini,fin))
        for tx_id, info in tx_dict.items():
            if info["strand"] == "+":
                info["exons"].sort(key=lambda x: x[0])
            else:
                info["exons"].sort(key=lambda x: x[0], reverse=True)
            long_tot = 0
            # El offset va a ser la longitud de exones acumulada antes del actual
            exons_offset = []
            for ini, fin in info["exons"]:
                long = fin - ini
                exons_offset.append((ini,fin,long_tot)) 
                long_tot += long
            info["exons_offset"]=exons_offset
            info["total_length"] = long_tot
    return tx_dict        

def genome_to_tx(ref_pos: int, tx_info: dict) -> int:
    '''
    We convert a genomic position to a transcript position based on the exons
    structure
    '''
    strand = tx_info["strand"]
    # Recorremos todos los exones
    for ini, fin, offset in tx_info["exons_offset"]:
        # Si la posicion está en el exon 
        if ini <= ref_pos < fin:
            # Si es pos 
            if strand == "+":
                return offset + (ref_pos - ini)
            else:
                # -1 porque la posicion fin no es inclusive
                return offset + (fin - 1 - ref_pos)
    return None
    
def tablaMods(args):
    ''' 
    We generate a tsv file with the CIGAR and Probs for each read in the bam 
    file. 
    '''
    mods_codigos = {nombre_a_codigo[m] for m in args.mods if m in nombre_a_codigo}
    not_mods = set()
    for filepath in [args.bam, args.bed, args.gtf, args.assoc]:
        if not os.path.exists(filepath):
            print(f"Error: File not found -> {filepath}")
            sys.exit(1)
    print(f"Loading dependencies for {os.path.basename(args.bam)} ...")
    conf_sites = cargaBed(args.bed)
    read_tx = cargaRead_Transcrito(args.assoc)
    tx_dict = cargaGTF(args.gtf)
    
    cpus = int(os.environ.get('SLURM_CPUS_PER_TASK', 4))
    print(f"Processing BAM utilizing {cpus} threads...")
    
    samfile = pysam.AlignmentFile(args.bam, "rb", threads=6)
    if not samfile.check_index():
        print("Error: The BAM file does not have an index (.bai). Aborting.")
        sys.exit(1)
    filas_lecturas = []
    
    for read in samfile.fetch():
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        
        read_id = read.query_name
        # Sacamos el transcrito al que está asociada la read
        tx_id = read_tx.get(read_id)
        if not tx_id or tx_id not in tx_dict:
            continue
        tx_info = tx_dict[tx_id]
        tx_length = tx_info["total_length"]
        
        length = read.query_length
        blocks = read.get_blocks()
        tx_start = genome_to_tx(blocks[0][0], tx_info)
        tx_end = genome_to_tx(blocks[-1][1] - 1, tx_info)
        
        if tx_start is None or tx_end is None:
            continue
        
        if tx_start > tx_end:
            tx_start, tx_end = tx_end, tx_start

        # Nuestra fila para el TSV final
        fila = {"isoform": read_id, 
                "transcript_id":tx_id,
                "transcript_length":tx_length,
                "tx_start":tx_start,
                "tx_end":tx_end,}
        
        try:
            # Dict[(canonical base, strand, modification)] -> [ (pos,qual), …] 
            mods = read.modified_bases
        except (ValueError, RuntimeError):
            continue
        last_pos = -1
        if mods:
            query_ref = dict(read.get_aligned_pairs(matches_only=True))
            chrom = read.reference_name
            strand = "-" if read.is_reverse else "+"
            valid_sites = conf_sites.get((chrom, strand), set())
            prob_threshold = args.prob_lim * 256
            # Para cada modificacion (canonical base, strand, modification) con su lista de [ (pos,qual), …] 
            
            for tupla, probs_list in mods.items():
                base_canon = tupla[0] 
                mod = str(tupla[2])
                if mod not in mods_codigos:
                    not_mods.add(mod)
                    continue
                # Para cuando el id de la mod es un codigo CHEBI
                mod_str = "M"
                valid_mods = []
                # Vamos a pasar las coordenadas de prob list a coords a nivel de transcrito
                for read_pos, prob in probs_list:
                    if prob < prob_threshold:
                        continue
                    
                    #Filtramos para asegurarnos de que está entre las posiciones de confianza
                    ref_pos = query_ref.get(read_pos)
                    if ref_pos is None or ref_pos not in valid_sites:
                        continue
                    # Sacamos la posicion a  nivel de transcrito
                    tx_pos = genome_to_tx(ref_pos,tx_info)
                    if tx_pos is not None:
                        valid_mods.append((tx_pos, prob))
                
                # We make sure that the list is sorted
                valid_mods.sort(key=lambda x: x[0])
                # Construimos el CIGAR
                cigar_parts = []
                prob_parts = []
                last_pos = -1
                
                for tx_pos, prob in valid_mods:
                    prob_frac = prob_ASCII(prob)
                    # Dist de unmodified (u)
                    dist_u = tx_pos -  last_pos -1
                    if dist_u > 0:
                        cigar_parts.append(f"{dist_u}U")
                    # Si la siguiente base mod va consecutiva
                    if dist_u == 0 and cigar_parts and cigar_parts[-1].endswith(mod_str):
                        count_previo = int(cigar_parts[-1][:-len(mod_str)])
                        cigar_parts[-1] = f"{count_previo + 1}{mod_str}"
                    else:
                        cigar_parts.append(f"1{mod_str}")
                    prob_parts.append(f"{prob_frac}")
                    last_pos = tx_pos
                bases_res = tx_length - last_pos -1
                
                if bases_res > 0:
                    cigar_parts.append(f"{bases_res}U")
                '''
                if prob_parts:
                    fila[f"CIGAR_{cod_mod_largo[mod]}"] = "".join(cigar_parts)
                    fila[f"Probs_{cod_mod_largo[mod]}"] = "".join(prob_parts)
                else:
                    fila[f"CIGAR_{cod_mod_largo[mod]}"] = f"{tx_length}U"
                    fila[f"Probs_{cod_mod_largo[mod]}"] = ""
                '''
                fila[f"positions_{cod_mod_largo[mod]}"] = [pos for pos, _ in valid_mods]
                fila[f"probabilities_{cod_mod_largo[mod]}"] = [prob_ASCII(prob) for _, prob in valid_mods]

            filas_lecturas.append(fila)
    samfile.close()
    # Guardar resultados
    df = pd.DataFrame(filas_lecturas)
    
    out_dir = os.path.dirname(args.out_tsv)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        
    df.to_csv(args.out_tsv, sep="\t", index=None)
    return f"Extraction completed. {len(df)} reads with modifications saved in {args.out_tsv}. Mods no cogidas: {sorted(not_mods)}" 


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Extrae tags de modificaciones a coordenadas de transcrito (1 muestra).")
    
    # Reemplazamos todos los directorios por rutas de archivos específicas
    parser.add_argument('--bam', required=True, help="Path to the BAM file")
    parser.add_argument('--bed', required=True, help="Path to the filtered BED file")
    parser.add_argument('--gtf', required=True, help="Path to the GTF annotation")
    parser.add_argument('--assoc', required=True, help="Path to the TSV of reads to isoforms from IsoQuant")
    parser.add_argument('--out_tsv', required=True, help="Path where the output TSV will be saved")
    parser.add_argument('--mods', nargs="+", required=True, help="Modifications to include")
    parser.add_argument('--prob_lim', type=float, default=0.95, help="Minimum probability threshold for modifications (default: 0.95)")
    
    args = parser.parse_args()
    
    print(tablaMods(args))