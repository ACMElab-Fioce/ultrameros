import os
import pandas as pd
import subprocess

# Define a pasta com os arquivos BAM
input_folder = "."  # Diretório atual onde o script está

# Lista para armazenar os resultados
results = []

# Loop para processar todos os arquivos .bam
for bam_file in os.listdir(input_folder):
    if bam_file.endswith("_ordenado.bam"):
        # Deriva o nome do FASTQ a partir do BAM
        base_name = bam_file.replace("_ordenado.bam", "")
        r1_file = f"{base_name}_R1_001.fastq.gz"
        r2_file = f"{base_name}_R2_001.fastq.gz"

        # Comando para obter o cabeçalho e extrair o nome da referência
        command_header = f"samtools view -H {bam_file}"
        try:
            header_output = subprocess.check_output(command_header, shell=True).decode('utf-8')
            # Extraindo o nome da referência
            references = []
            for line in header_output.splitlines():
                if line.startswith("@SQ"):
                    ref_info = line.split('\t')
                    ref_name = [info for info in ref_info if info.startswith("SN:")]
                    if ref_name:
                        references.append(ref_name[0].split(":")[1])
            
            # Para cada referência, calcular o número de reads alinhadas
            for ref in references:
                command_idxstats = f"samtools idxstats {bam_file} | grep {ref}"
                output = subprocess.check_output(command_idxstats, shell=True).decode('utf-8')
                aligned_reads_count = sum(int(line.split()[2]) for line in output.strip().split('\n') if line)
                
                # Adiciona os resultados à lista com o nome do FASTQ
                results.append({'Reference': ref, 'Aligned Reads': aligned_reads_count, 'FASTQ': r1_file})

        except subprocess.CalledProcessError as e:
            print(f"Erro ao processar {bam_file}: {e}")

# Cria um DataFrame e salva como CSV
df = pd.DataFrame(results)
df.to_csv("aligned_reads_summary.csv", sep=';', index=False)

print("Resumo salvo em aligned_reads_summary.csv")
