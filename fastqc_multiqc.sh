#!/bin/bash

#Este script permite realizar los controles de calidad para las muestras preprocesadas o procesadas (se adaptó según el procesamiento se había realizado o no)

#-----------------------------------------------------------------------------------------
# ANÁLISIS DE CALIDAD DE LAS MUESTRAS CON FASTQC Y MULTIQC
#-----------------------------------------------------------------------------------------

module load fastqc


# Análisis de calidad por muestra + total de muestras 

# Análisis de todos los archivos fastqs:

for file in *.fastq; do
    if [[ -f "$file" ]]; then  
        echo "Procesando $file con FastQC..."
        fastqc -o fastqc_reports "$file"  
    fi
done


# MultiQC

echo "Generando informe con MultiQC..."
multiqc fastqc_reports -o multiqc_report

echo "Análisis completado. Los informes están en las carpetas fastqc_reports y multiqc_report"
