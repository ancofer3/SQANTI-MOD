import pysam
import pandas as pd
import glob
import os
import gzip
from concurrent.futures import ProcessPoolExecutor


# A futuro esto deberíamos cambiarlo desde el principio o dar tu un dic
PROB_LIM=0.5
cod_mod_largo = {"m":"m5C","a":"m6A", "19228":"2OmeC", "17802":"pseU","19227":"2OmeU", "17596":"inosine" , "69426":"2OmeA"}
cod_mod_1 = {"m":"m","a":"a", "19228":"C", "17802":"P","19227":"U", "17596":"I" , "69426":"A"}

def prob_ASCII(prob:float) -> str:
    if prob >= 0 and prob <= 256:
        ascii_code = chr(33 + round((prob * 93) / 256))
    else:
        print(prob)
        raise ValueError
    return ascii_code

def cargaBed(path: str) -> set:
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

def cargaRead_Transcrito(path:str) -> dict:
    # Sacamos la relacion read_id-transcript_id
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
    "Parseamos GTF para sacar la estructura de exones y al longitud del transcrito"
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
    
def tablaMods(path):
    nombre = "_".join(os.path.basename(path).split(".")[0].split("_")[0:2])
    
    # Cargamos toda la información complementaria a los sams
    bed_path = f"data/beds/{nombre}_filtered.bed" 
    gtf_path = f"reconstruction/isoquant/{nombre}/{nombre}/{nombre}.extended_annotation.gtf"
    assoc_path = f"reconstruction/isoquant/{nombre}/{nombre}/{nombre}.transcript_model_reads.tsv.gz"
    if not os.path.exists(gtf_path) or not os.path.exists(assoc_path):
        return f"Error: Archivos de IsoQuant no encontrados para {nombre}"
    
    conf_sites = cargaBed(bed_path)
    read_tx = cargaRead_Transcrito(assoc_path)
    tx_dict = cargaGTF(gtf_path)
    print("El set conf_sites tiene esta pinta:", list(conf_sites)[0:5])
    
    samfile = pysam.AlignmentFile(path, "rb", threads=6)
    if not samfile.check_index():
        raise ValueError
    filas_lecturas = []
    
    for read in samfile.fetch():
        if read.is_secondary or read.is_supplementary or read.query_sequence is None:
            continue
        
        read_id = read.query_name
        # Sacamos el transcrito al que está asociada la read
        tx_id = read_tx.get(read_id)
        if not tx_id or tx_id not in tx_dict:
            continue
        tx_info = tx_dict[tx_id]
        tx_length = tx_info["total_length"]
        
        length = read.query_length
        # Nuestra fila para el TSV final
        fila = {"isoform": read_id, 
                "transcript_id":tx_id,
                "transcript_length":tx_length}
        
        try:
            # Dict[(canonical base, strand, modification)] -> [ (pos,qual), …] 
            mods = read.modified_bases
        except (ValueError, RuntimeError):
            continue
        last_pos = -1
        if mods:
            # Hacemos un dic de posicion en la read a genomica
            query_ref = {q:r for q,r in read.get_aligned_pairs(matches_only=True) if q is not None and r is not None}
            chrom = read.reference_name
            strand = "-" if read.is_reverse else "+"
            # Para cada modificacion (canonical base, strand, modification) con su lista de [ (pos,qual), …] 
            for tupla, probs_list in mods.items():
                base_canon = tupla[0] 
                mod = str(tupla[2])
                # Para cuando el id de la mod es un codigo CHEBI
                mod_str = "m"
                valid_mods = []
                # Vamos a pasar las coordenadas de prob list a coords a nivel de transcrito
                for read_pos, prob in probs_list:
                    if (prob/256) < PROB_LIM:
                        continue
                    
                    #Filtramos para asegurarnos de que está entre las posiciones de confianza
                    ref_pos = query_ref.get(read_pos)
                    if ref_pos is None or (chrom, ref_pos, strand) not in conf_sites:
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
                if prob_parts:
                    fila[f"CIGAR_{cod_mod_largo[mod]}"] = "".join(cigar_parts)
                    fila[f"Probs_{cod_mod_largo[mod]}"] = "".join(prob_parts)
                else:
                    fila[f"CIGAR_{cod_mod_largo[mod]}"] = f"{tx_length}U"
                    fila[f"Probs_{cod_mod_largo[mod]}"] = ""

            filas_lecturas.append(fila)
    samfile.close()
    df = pd.DataFrame(filas_lecturas)
    os.makedirs("tsvs_transcritos", exist_ok=True)
    df.to_csv(f"tsvs_transcritos/{nombre}_transcripts_modCIGAR.tsv",sep="\t",index=None)
    return f"{nombre} Completado con éxito. {len(df)} lecturas procesadas" 


if __name__ == '__main__':
    cpus_slurm = int(os.environ.get('SLURM_CPUS_PER_TASK', 1))
    archivos_bam = glob.glob("data/bams_primary/*primary.bam")
    print(f"Se han encontrado {len(archivos_bam)} archivos BAM.")
    print(f"Ejecutando en PARALELO utilizando {cpus_slurm} CPUs...\n")
    # Lanzar un "Pool" de trabajadores. Procesará N archivos a la vez.
    with ProcessPoolExecutor(max_workers=len(archivos_bam)) as executor:
        resultados = executor.map(tablaMods, archivos_bam)
        # Imprime los avisos de finalización conforme van acabando
        for res in resultados:
            print(res)
    print("\n¡Todos los archivos han sido procesados!")
