#!/bin/bash

#Este script permite la cuantificación con HTSeq y procesamiento de lecturas

#Para el conteo de exones se utilizó el comando ‘htseq-count’ en modo ‘intersection-strict’
#Para el conteo total, se empleó el modo ‘union’ para permitir solapamientos parciales con ciertas estructuras, como los intrones.


#-------------------------------------------------------------------------------------------
#UNION
#------------------------------------------------------------------------------------------

#Nota: donde está el entrecomillado se utilizaron las rutas a cada uno de los directorios (el usuario debe meter las suyas propias)
bam_dir=""
gtf_file=""
output_dir=""


muestras=($bam_dir/SRR220022{65..99}_Aligned.out_sorted.bam $bam_dir/SRR220023{00..30}_Aligned.out_sorted.bam)

bam_file="${muestras[$SLURM_ARRAY_TASK_ID]}"

# Recorrer todos lois archivos .bam en el directorio bam_dir
if [[ -f "$bam_file" ]]; then
    # Obtener el nombre base del archivo (sin ruta ni extensión)
    base_name=$(basename "$bam_file" .bam)

    # Definir el archivo de salida
    output_file="$output_dir/${base_name}_htseq.txt"

    # Ejecutar htseq-count para cada archivo BAM
    echo "Procesando $bam_file ..."
    htseq-count -f bam -m union "$bam_file" "$gtf_file" > "$output_file"

    echo "Resultado guardado en $output_file"
fi



#--------------------------------------------------------------------------
#INTERSECTION-STRICT
#---------------------------------------------------------------------------

bam_dir=""
gtf_file=""
output_dir=""

muestras=($bam_dir/SRR22002{247..330}_Aligned.out_sorted.bam)
bam_file="${muestras[$SLURM_ARRAY_TASK_ID]}"

# Recorrer todos los archivos .bam en el directorio bam_dir


if [[ -f "$bam_file" ]]; then
    # Obtener el nombre base del archivo (sin ruta ni extensión)
    base_name=$(basename "$bam_file" .bam)

    # Definir el archivo de salida
    output_file="$output_dir/${base_name}_htseq.txt"

    # Ejecutar htseq-count para cada archivo BAM
    echo "Procesando $bam_file ..."
    htseq-count -f bam -m intersection-strict "$bam_file" "$gtf_file" > "$output_file"

    echo "Resultado guardado en $output_file"
fi



#---------------------------------------------------------------------------
#CÁLCULO DE LAS LECTURAS INTRÓNICAS
#---------------------------------------------------------------------------

dir_total=".../output_htseq/total_reads"
dir_exons=".../output_htseq/exon_reads"
dir_introns=".../output_htseq/introns_reads"

# Crear directorio de salida si no existe
mkdir -p "$dir_introns"

# Lista de muestras
muestras=()
muestras=("SRR22002299")


# Seleccionar archivo basado en el ID de la tarea
sample_id="${muestras[$SLURM_ARRAY_TASK_ID]}"

file_total="$dir_total/${sample_id}_Aligned.out_sorted_htseq.txt"
file_exons="$dir_exons/${sample_id}_Aligned.out_sorted_htseq.txt"
file_introns="$dir_introns/${sample_id}_introns.txt"

# Verificar si los archivos existen

if [[ -f "$file_total" && -f "$file_exons" ]]; then
    echo "Procesando muestra: $sample_id"

    # Calcular lecturas intrónicas
    awk '
        NR==FNR { total[$1] = $2; next }
        $1 in total {
            if ($1 ~ /^__/) next; # Ignorar líneas especiales
            intron = total[$1] - $2;
            print $1, intron;
        }
    ' "$file_total" "$file_exons" > "$file_introns"

    echo "Archivo procesado: $file_introns"
else
    echo "Advertencia: No se encontró el archivo total o exons para '$sample_id'."
    exit 1
fi
