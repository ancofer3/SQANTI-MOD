import pandas as pd
import glob
import re

# Aqui tenemos que hacer que saque una occupancy local en función de cuantas reads solapan con esa posicion

def collapseCIGARs(tsv):
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
		# Para cada modificacion posible
		for mod in mod_types:
			positions_col = f"positions_{mod}"
			pos_counts = {}
			for positions in group[positions_col].dropna():
				for pos in positions:
					if pos not in pos_counts:
						pos_counts[pos] = 0
					pos_counts[pos] += 1
			# Nos quedamos con todas las posiciones posibles
			sorted_pos = sorted(pos_counts.keys())
			cigar_parts = []
			occupancies = []
			
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
