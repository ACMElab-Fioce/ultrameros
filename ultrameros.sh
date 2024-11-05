#!/bin/bash
## Job name
#PBS -N LoteUltrameros

# Redirect output stream to this file.
#PBS -o pbs_out.dat

# Redirect error stream to this file.
#PBS -e pbs_err.dat

### Number of nodes and number of processors per node
#PBS -l nodes=1:ppn=96

# Change to current working directory (directory where qsub was executed)
cd $PBS_O_WORKDIR

# Conda environment setup
__conda_setup="$('/home/pjeronimo/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/pjeronimo/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/home/pjeronimo/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/pjeronimo/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup

# Ativar ambiente conda
conda activate teste-ultrameros

# Caminho do arquivo de referência
reference="Article_core_seq.fasta"

# Função para processar cada par de arquivos FASTQ
processar_bam() {
    local r1="$1"
    local r2="$2"
    local base="$3"

    # Controle de qualidade usando fastp
    fastp -i "$r1" -I "$r2" -o "${base}_cleaned_R1.fastq.gz" -O "${base}_cleaned_R2.fastq.gz" -q 20

    # Alinhar com BWA e gerar arquivo SAM
    bwa mem $reference "${base}_cleaned_R1.fastq.gz" "${base}_cleaned_R2.fastq.gz" > ${base}.sam

    # Converter SAM para BAM
    samtools view -S -b ${base}.sam > ${base}.bam

    # Ordenar o arquivo BAM
    samtools sort -o ${base}_ordenado.bam ${base}.bam

    # Indexar o arquivo BAM
    samtools index ${base}_ordenado.bam
}

# Loop para processar todos os arquivos fastq.gz
for r1 in *R1_001.fastq.gz; do
    r2="${r1/_R1_001/_R2_001}"
    base="${r1%%_R1_001.fastq.gz}"

    # Chamar a função para processar o par R1 e R2
    processar_bam $r1 $r2 $base
done

echo "Processamento concluído."
