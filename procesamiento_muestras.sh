#!/bin/bash

#Este script parte de los datos descargados del repositorio GEO.
#Se descargó el archivo de texto 'Accession List', a partir del cual se obtuvieron los archivos '.sra' para cada muestra y se extrajeron los fastqs de cada muestra


#-----------------------------------------------------------
#PROCESAMIENTO DE ARCHIVOS FASTQ
#-----------------------------------------------------------
#Procesamiento de fastqs: eliminación de adaptadores y de cola poli-A 
#eliminación de bases iniciales y procesado por calidad < 30 


#Directorios Input y output
input_dir=".../analysis/fastqs"
output_dir=".../analysis/fastqs_processed"
report_dir=".../analysis/reports_cutadapt"

mkdir -p "$output_dir"
mkdir -p "$report_dir"

#Iteración sobre los paired FASTQ en el directorio input

for file1 in "$input_dir"/*_1.fastq; do
        #Basename of file w/o suffixes
        base_name=$(basename "$file1" _1.fastq)
        file2="$input_dir/${base_name}_2.fastq"
        if [[ -f "$file2" ]]; then
                echo "Procesando archivos: $file1 y $file2"


                # Eliminación de adaptadores
                 cutadapt --cores=10 \
                        -a AGATCGGAAGAG \
                        -A AGATCGGAAGAG \
                        --minimum-length 20 \
                        -o "$output_dir/${base_name}_1_trim.fastq" \
                        -p "$output_dir/${base_name}_2_trim.fastq" \
                        "$file1" "$file2" \
                        > "$report_dir/report_${base_name}_trim.txt"



                # Eliminación de cola Poli-A
                cutadapt --cores=10 \
                        -a "AAAAAAAAAA" \
                        -A "AAAAAAAAAA" \
                        --minimum-length 20 \
                        -o "$output_dir/${base_name}_1_trim_polyA.fastq" \
                        -p "$output_dir/${base_name}_2_trim_polyA.fastq" \
                         "$output_dir/${base_name}_1_trim.fastq" "$output_dir/${base_name}_2_trim.fastq" \
                        > "$report_dir/report_${base_name}_trim_polyA.txt"

                # Eliminación de bases iniciales 
                cutadapt --cores=10 \
                        -u 25 \
                        -U 25 \
                        --minimum-length 20 \
                        -o "$output_dir/${base_name}_1_trim_polyA_trimmed.fastq" \
                        -p "$output_dir/${base_name}_2_trim_polyA_trimmed.fastq" \
                        "$output_dir/${base_name}_1_trim_polyA.fastq" "$output_dir/${base_name}_2_trim_polyA.fastq" \
                        > "$report_dir/report_${base_name}_trim_polyA_trimmed.txt"


                # Procesado por calidad
                cutadapt --cores=10 \
                        -q 30 -Q 30 -u 3 -U 3 \
                          --minimum-length 20 \
                        -o "$output_dir/${base_name}_1_trim_polyA_trimmed_Q.fastq" \
                        -p "$output_dir/${base_name}_2_trim_polyA_trimmed_Q.fastq" \
                        "$output_dir/${base_name}_1_trim_polyA_trimmed.fastq" "$output_dir/${base_name}_2_trim_polyA_trimmed.fastq" \
                        > "$report_dir/report_${base_name}_trim_polyA_trimmed_Q.txt"

        else
                echo "Archivo emparejado para $file1 no encontrado. Saltando..."
        fi
done

echo "Procesamiento completado"
