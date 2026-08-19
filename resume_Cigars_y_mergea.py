import pandas as pd
import glob
import re

def collapseCIGARs(tsv):
	cigar_cols = [col for col in tsv.columns if col.startswith("CIGAR")]
	mod_types = [col.replace("CIGAR_","") for col in cigar_cols]
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
		# Para cada modificacion posible
		for mod in mod_types:
			cigar_col = f"CIGAR_{mod}"
			pos_counts = {}
			for cigar in group[cigar_col].dropna():
				if cigar == f"{tx_length}U":
					continue
				tokens = re.findall(r'(\d+)([A-Za-z]+)', str(cigar))
				current_pos = 0
				# Para cada token (num,letra) del CIGAR
				for len_str, op in tokens:
					length = int(len_str)
					if op == "U":
						current_pos += length
					else:
						for i in range(length):
							real_pos = current_pos + i
							pos_counts[real_pos] = pos_counts.get(real_pos, 0) + 1
						current_pos += length
			# Nos quedamos con todas las posiciones posibles
			sorted_pos = sorted(pos_counts.keys())
			cigar_parts = []
			occupancies = []
			last_pos = -1
			# Remontamos el CIGAR
			for pos in sorted_pos:
				dist_u = pos -last_pos -1
				if dist_u > 0:
					cigar_parts.append(f"{dist_u}U")
				count = pos_counts[pos]
				occupancy_abs = count
				occupancies.append(occupancy_abs)
				if dist_u == 0 and cigar_parts and cigar_parts[-1].endswith("m"):
					count_previo = int(cigar_parts[-1][:-1])
					cigar_parts[-1] = f"{count_previo + 1}m"
				else:
					cigar_parts.append(f"1m")
				last_pos = pos
			bases_res = tx_length - last_pos - 1
			if bases_res > 0:
				cigar_parts.append(f"{bases_res}U")
			# Y lo añadimos
			fila_tx[f"CIGAR_{mod}"] = "".join(cigar_parts)
			fila_tx[f"Occupancy_{mod}"] = occupancies
		resultados_transcritos.append(fila_tx)
	df_transcrito = pd.DataFrame(resultados_transcritos)
	return df_transcrito

for i in glob.glob("SQANTI3_QC_isoquant_*/*classification.txt"):
	name = i.split("/")[-1].replace("_classification.txt","")
	print("Empezando a procesar:",name)
	cls = pd.read_csv(i,sep="\t")
	tsv = pd.read_csv(f"tsvs_transcritos/{name}_transcripts_modCIGAR.tsv",sep="\t")
	tsv_colapsado = collapseCIGARs(tsv)
	fusion = cls.merge(tsv_colapsado,how="left",left_on="isoform",right_on="transcript_id")
	fusion.to_csv(f"data/classifications_mod/{name}_classification_transcripts_mod.tsv", sep="\t")
	print("Mergeo terminado con: ",name)
