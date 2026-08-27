import ast
import pandas as pd
import numpy as np
import os
import argparse
import sys

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
			fila_tx[f"Positions_{mod}"] = ";".join(sorted_pos)
			fila_tx[f"Coverage_{mod}"] = ";".join(local_coverages)
			fila_tx[f"Occupancy_{mod}"] = ";".join(occupancies)
            
		resultados_transcritos.append(fila_tx)
	df_transcrito = pd.DataFrame(resultados_transcritos)
	return df_transcrito

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--tsv_in', required=True, help="Ruta al TSV con modificaciones de la muestra")
    parser.add_argument('--sqanti_class', required=True, help="Ruta al classification.txt de SQANTI3 para la muestra")
    parser.add_argument('--out_tsv', required=True, help="Ruta al archivo TSV de salida final")
    
    parser.add_argument('--min_cov', type=int, default=20, help="Cobertura mínima por posición")
    parser.add_argument('--min_occ', type=float, default=0.2, help="Ocupancia mínima por posición")
    args = parser.parse_args()
    
    print(f"Empezando a procesar archivo: {args.tsv_in}")
    
    # Comprobar que los archivos de entrada existen
    if not os.path.exists(args.tsv_in):
        print(f"Error: input {args.tsv_in} not found")
        sys.exit(1)
        
    if not os.path.exists(args.sqanti_class):
        print(f"Error: SQANTI classification file not found: {args.sqanti_class}")
        sys.exit(1)
        
    # Crear carpeta de destino si no existe
    out_dir = os.path.dirname(args.out_tsv)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        
    # 1. Leer y colapsar el TSV de modificaciones
    tsv = pd.read_csv(args.tsv_in, sep="\t")
    tsv_colapsado = collapsePositions(tsv, args)
    
    # 2. Leer clasificación de SQANTI
    cls = pd.read_csv(args.sqanti_class, sep="\t")
    
    # 3. Mergear
    fusion = cls.merge(tsv_colapsado, how="left", left_on="isoform", right_on="transcript_id")
    
    # 4. Guardar archivo final
    fusion.to_csv(f"{args.out_tsv}", sep="\t", index=False)
    
    print(f"Merging finished. Results stored in {out_dir}/{args.out_tsv}")

