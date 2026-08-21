import ast
import pandas as pd
import glob
import numpy as np
import os
import argparse

'''OJO: Faltaria ver como incorporamos las probabilidades. Alomejor lo mejor es que ya desde el principio se de
la prob límite que quieres usar para cada mod y en base a eso extraigamos antes solo las que pasen el filtro.
'''
# FILTERS recommended by Ghohabi Esfahani et al. (2026) (DOI: 10.1186/s13059-026-04096-w)
FILTRO_cov = 20
FILTRO_occ = 0.2
def pos_to_CIGAR(positions, tx_length):
    # Retorno rápido si no hay modificaciones
    if not positions:
        return f"{tx_length}U"

    # Eliminamos posibles duplicados y ordenamos
    sorted_positions = sorted(list(set(positions)))
    cigar_parts = []
    
    # Puntero para rastrear el final de la última anotación (asumiendo coordenadas 0-based)
    last_unmodified_end = 0 
    
    current_start = sorted_positions[0]
    current_length = 1
    
    # Agrupación de posiciones
    for i in range(1, len(sorted_positions)):
        if sorted_positions[i] == sorted_positions[i - 1] + 1:
            current_length += 1
        else:
            # 1. Rellenar el gap previo de 'U' si el bloque 'M' no empieza justo donde acabó el anterior
            if current_start > last_unmodified_end:
                cigar_parts.append(f"{current_start - last_unmodified_end}U")
            
            # 2. Insertar el bloque modificado 'M'
            cigar_parts.append(f"{current_length}M")
            
            # 3. Actualizar el puntero y resetear contadores
            last_unmodified_end = current_start + current_length
            current_start = sorted_positions[i]
            current_length = 1

    # Procesar el último bloque 'M' remanente fuera del loop
    if current_start > last_unmodified_end:
        cigar_parts.append(f"{current_start - last_unmodified_end}U")
    
    cigar_parts.append(f"{current_length}M")
    
    # Calcular y añadir el sufijo 'U' si no hemos alcanzado el límite del transcrito
    last_unmodified_end = current_start + current_length
    if last_unmodified_end < tx_length:
        cigar_parts.append(f"{tx_length - last_unmodified_end}U")
        
    return "".join(cigar_parts)
    
def collapsePositions(tsv, args):
	pos_cols = [col for col in tsv.columns if col.startswith("positions_")]
	mod_types = [col.replace("positions_","") for col in pos_cols]
	resultados_transcritos = []
	# Para cada transcript_id
	for tx_id, group in tsv.groupby("transcript_id"):
		tx_length = group["transcript_length"].iloc[0]
		total_reads = len(group)
		fila_tx = {
			"transcript_id": tx_id,
			"transcript_length": tx_length,
			"total_reads" : total_reads
		}
		# Asegurar que las coordenadas son numéricas
		starts = group["tx_start"].astype(int).values
		ends = group["tx_end"].astype(int).values
		# Para cada modificacion posible
		for mod in mod_types:
			positions_col = f"positions_{mod}"
			pos_counts = {}
			for positions in group[positions_col].dropna():
				if isinstance(positions, str):
					try:
						if positions.startswith('['):
							positions = ast.literal_eval(positions)
						else:
							positions = positions.split(',')
					except (ValueError, SyntaxError):
						continue
				for pos in positions:
					pos = int(pos)
					if pos in pos_counts:
						pos_counts[pos] += 1
					else:
						pos_counts[pos] = 1
			# Vamos a quitarnos todas las posiciones que no tengan suficiente cobertura
			pos_filt = {pos: count for pos, count in pos_counts.items() if count >= args.min_cov}
			# Tambien nos quitamos las que no cumplan con el filtro de ocupancia
			pos_filt = {pos: count for pos, count in pos_filt.items() if count / np.sum((starts <= pos) & (ends >= pos)) >= args.min_occ}
			sorted_pos = sorted(pos_filt.keys())
			
			local_coverages = []
			occupancies = []

			for pos in sorted_pos:
				# Una lectura cubre pos si tx_start <= pos <= tx_end
				n_cov = int(np.sum((starts <= pos) & (ends >= pos)))
				n_mod = pos_counts[pos]
				occ = float(round(n_mod / n_cov, 4) if n_cov > 0 else 0.0)

				local_coverages.append(n_cov)
				occupancies.append(occ)
			
			# Y lo añadimos
			fila_tx[f"Positions_{mod}"] = sorted_pos
			fila_tx[f"Coverage_{mod}"] = local_coverages
			fila_tx[f"Occupancy_{mod}"] = occupancies
		resultados_transcritos.append(fila_tx)
	df_transcrito = pd.DataFrame(resultados_transcritos)
	return df_transcrito

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--tsv_dir', required=True, help="Directorio con los TSV de entrada")
    parser.add_argument('--sqanti_dir', required=True, help="Directorio padre con los outputs de SQANTI3")
    parser.add_argument('--out_dir', required=True, help="Directorio de salida final")
    parser.add_argument('--min_cov', type=int, default=20, help="Cobertura mínima por posición")
    parser.add_argument('--min_occ', type=float, default=0.2, help="Ocupancia mínima por posición")
    args = parser.parse_args()
	
    os.makedirs(args.out_dir, exist_ok=True)
    patron_sqanti = os.path.join(args.sqanti_dir, "SQANTI3_QC_isoquant_*", "*classification.txt")
    for i in glob.glob(patron_sqanti):
        name = os.path.basename(i).replace("_classification.txt", "")
        print("Empezando a procesar:", name)
        
        cls = pd.read_csv(i, sep="\t")
        tsv_path = os.path.join(args.tsv_dir, f"{name}_transcripts_modCIGAR.tsv")
        
        if not os.path.exists(tsv_path):
            print(f"Aviso: No se encontró {tsv_path}, saltando...")
            continue
            
        tsv = pd.read_csv(tsv_path, sep="\t")
        tsv_colapsado = collapsePositions(tsv, args)
        
        fusion = cls.merge(tsv_colapsado, how="left", left_on="isoform", right_on="transcript_id")
        out_file = os.path.join(args.out_dir, f"{name}_classification_transcripts_mod.tsv")
        fusion.to_csv(out_file, sep="\t", index=False)
        
        print("Mergeo terminado con: ", name)

