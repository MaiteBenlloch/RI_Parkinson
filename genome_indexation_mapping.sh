#!/bin/bash

#Este script permite la indexación y el mapeo contra el genoma de referencia una vez se descargó.
#Además, también permite ordenar los archivos BAM e indexarlos
#Partimos de los archivos de genoma y anotación:
#Homo_sapiens.GRCh38.113.gtf
#Homo_sapiens.GRCh38.dna.primary_assembly.fa

#----------------------------------------------------------------------------------
#INDEXACIÓN
#----------------------------------------------------------------------------------
STAR --runMode genomeGenerate --runThreadN 2 --genomeDir genome --genomeFastaFiles genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile genome/Homo_sapiens.GRCh38.113.gtf --sjdbOverhang 149 



#---------------------------------------------------------------------------------
#MAPEO CONTRA GENOMA DE REFERENCIA
#---------------------------------------------------------------------------------
# directorio del genoma
genome_dir="genome/"

# directorio de las muestras fastq
input_dir="fastqs_processed/"


output_dir="output_STAR/"

# crear el directorio de salida si no existe
mkdir -p "$output_dir"

# iterar sobre los pares de reads fastq (_1 y _2)
for r1 in "$input_dir"/*_1_trim_polyA_trimmed_Q.fastq; do
    # obtener el nombre base del archivo r1
    basename=$(basename "$r1" _1_trim_polyA_trimmed_Q.fastq)

    # archivo r2 correspondiente
    r2="$input_dir/${basename}_2_trim_polyA_trimmed_Q.fastq"

    # verificar si el archivo r2 existe
    if [[ -f "$r2" ]]; then
        echo "Procesando muestra: $basename"

        # ejecutar star
        STAR --runThreadN 4 \
             --genomeDir "$genome_dir" \
             --readFilesIn "$r1" "$r2" \
             --quantMode GeneCounts \
             --outFileNamePrefix "$output_dir/${basename}_" \
             --outSAMtype BAM Unsorted
    else
        echo "Archivo R2 no encontrado para $basename. Omitiendo."
    fi
done

echo "Mapeo completado."


#---------------------------------------------------------------------------
#ORDENAMIENTO E INDEXACIÓN DE LOS ARCHIVOS BAM
#---------------------------------------------------------------------------

input_dir="output_STAR/"
output_dir="sorted_bam/"

mkdir -p "$output_dir"

# iterar sobre los archivos bam en el directorio de entrada
for bam_file in "$input_dir"/*.bam; do
    # obtener el nombre base del archivo bam
    basename=$(basename "$bam_file" .bam)

    echo "Procesando $basename..."

    # archivo de salida ordenado
    sorted_bam="$output_dir/${basename}_sorted.bam"

    # ordenar el archivo bam
    samtools sort -o "$sorted_bam" "$bam_file"

    # indexar el archivo bam ordenado
    samtools index "$sorted_bam"

    echo "Archivo ordenado e indexado: $sorted_bam"
done
echo "Procesamiento completado."

