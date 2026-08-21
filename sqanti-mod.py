import argparse
import os
import subprocess
import sys
import glob

def run_command(comando, descripcion):
    """Función auxiliar para ejecutar comandos en la terminal desde Python"""
    print(f"\n[SQANTI-MOD] Iniciando: {descripcion}")
    print(f"Comando: {' '.join(comando)}")
    
    # Ejecuta el comando y detiene el pipeline si hay un error
    proceso = subprocess.run(comando)
    if proceso.returncode != 0:
        print(f"Fatal error executing {descripcion}. Aborting pipeline.")
        sys.exit(1)
    print(f"{descripcion} completed succesfully.\n")

def main():
    parser = argparse.ArgumentParser(description="SQANTI-MOD V1.0: Módulo independiente para modificaciones de ARN")
    
    # Inputs principales (Los que tú has definido)
    parser.add_argument('--bams', required=True, help="Carpeta con los BAMs modificados")
    parser.add_argument('--sqanti_class', required=True, help="Carpeta con los classifications de SQANTI3")
    parser.add_argument('--gtfs', required=True, help="Carpeta con los GTFs de reconstrucción (IsoQuant)")
    parser.add_argument('--read_assocs', required=True, help="Carpeta con los TSVs que relacionan reads e isoformas")
    parser.add_argument('--out_dir', required=True, help="Carpeta de salida final")
    # Parámetros para Modkit
    parser.add_argument('--run_modkit', action='store_true', help="Si se activa, el pipeline ejecuta Modkit internamente")
    parser.add_argument('--modkit_bin', default="modkit", help="Ruta al ejecutable de modkit")
    parser.add_argument('--ref_genome', help="Genoma de referencia (Obligatorio si usas --run_modkit)")
    
    args = parser.parse_args()

    # 0. Crear directorios temporales y de salida
    beds_dir = os.path.join(args.out_dir, "01_beds")
    tsvs_dir = os.path.join(args.out_dir, "02_extracted_mods")
    final_dir = os.path.join(args.out_dir, "03_final_classification")
    
    for d in [beds_dir, tsvs_dir, final_dir]:
        os.makedirs(d, exist_ok=True)

    # 1. PASO OPCIONAL: Ejecutar Modkit
    # (Si el usuario ya tiene los BEDs, se saltaría este paso)
    if args.run_modkit:
        if not args.ref_genome:
            print("Error: Para correr modkit necesitas proporcionar --ref_genome")
            sys.exit(1)
            
        bams = glob.glob(os.path.join(args.bams, "*.bam"))
        for bam in bams:
            nombre = os.path.basename(bam).replace(".bam", "")
            cmd_modkit = [
                args.modkit_bin, "pileup", bam, os.path.join(beds_dir, f"{nombre}.bed"),
                "--reference", args.ref_genome,
                "--modified-bases", "m5C", "m6A" # Puedes parametrizar esto también
            ]
            run_command(cmd_modkit, f"Modkit pileup para {nombre}")
            
            # Aquí podrías añadir el comando awk de filtrado con subprocess
    else:
        print("[SQANTI-MOD] Saltando Modkit (asumiendo que los BEDs ya existen o se proveen)")
        # En una versión más avanzada, pedirías un '--beds_dir' si no corres modkit.

    # 2. PASO: Extraer mods a coordenadas de transcrito
    # Llamamos a tu script refactorizado
    cmd_extract = [
        "python", "saca_tags_transcript_models.py",
        "--bam_dir", args.bams,
        "--bed_dir", beds_dir, 
        "--isoquant_dir", args.gtfs, # Adaptar según la estructura final que decidas
        "--out_dir", tsvs_dir
    ]
    run_command(cmd_extract, "Extracción de modificaciones a transcritos")

    # 3. PASO: Colapsar y mergear con SQANTI
    cmd_merge = [
        "python", "collapse_to_tx_and_merge.py",
        "--tsv_dir", tsvs_dir,
        "--sqanti_dir", args.sqanti_class,
        "--out_dir", final_dir
    ]
    run_command(cmd_merge, "Fusión final con SQANTI3 classification")

    print(f"Pipeline SQANTI-MOD completado. Resultados en: {final_dir}")

if __name__ == '__main__':
    main()