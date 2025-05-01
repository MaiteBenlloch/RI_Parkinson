# Análisis de Retención Intrónica en la Enfermedad de Parkinson

Este repositorio contiene los scripts y el flujo de trabajo utilizados para analizar el porcentaje de **retención intrónica (IR)** en muestras de pacientes con enfermedad de Parkinson, clasificadas según las **etapas de Braak**, que describen la progresión neuropatológica de la enfermedad.

## Descripción general del análisis

El objetivo es calcular la retención intrónica por muestra y por etapa de Braak, identificar diferencias significativas entre grupos y realizar un análisis de enriquecimiento funcional de los genes afectados.

---

##  Flujo de trabajo

### 1. Análisis de calidad inicial

- **Script:** `fastqc_multiqc.sh`
- **Descripción:** Control de calidad de los archivos FASTQ utilizando FastQC y MultiQC.
- **Nota:** Puede adaptarse tanto a muestras crudas (preprocesadas) como a muestras ya procesadas, según el diseño experimental.

### 2. Preprocesamiento de muestras

- **Script:** `procesamiento_muestras.sh`
- **Descripción:**
  - Eliminación de adaptadores
  - Eliminación de la cola poli-A
  - Recorte de bases iniciales
  - Quality trimming

### 3. Análisis de calidad posterior al procesamiento

- **Script:** `fastqc_multiqc.sh`
- **Descripción:** Se vuelve a ejecutar FastQC y MultiQC para confirmar la mejora en la calidad tras el preprocesamiento.

### 4. Descarga del genoma, mapeo y ordenamiento e indexación de archivos 

- **Script:** `genome_indexation_mapping.sh`
- **Descripción:**
  - Descarga e indexación del genoma de referencia (con STAR)
  - Mapeo de las lecturas (con STAR)
  - Ordenación e indexación de los archivos BAM
Nota: es importante ajustar el parámetro --sjdbOverhang a los datos del usuario (longitud de lectura - 1)

### 5. Cuantificación de intrones

- **Script:** `HTSeq.sh`
- **Descripción:** Cuantificación de lecturas que mapean a intrones usando HTSeq.

### 6. Cálculo del porcentaje de retención intrónica y análisis en R

- **Script:** `PercIR_GeneAnalysis.R`
- **Descripción:**
  - Cálculo del porcentaje de retención intrónica por muestra y etapa de Braak
  - Análisis estadístico de diferencias significativas (Kruskal-Wallis)
  - Visualización con boxplots del IR% por etapa
  - Análisis de enriquecimiento de rutas biológicas para genes con alta retención intrónica

---

##  Requisitos

- **Software:**
  - R (y RStudio)
  - FastQC, MultiQC
  - Cutadapt
  - STAR
  - HTSeq
  - Este trabajo ha sido realizado mediante conexión remota desde un terminal WSL a un clúster de computación de alto rendimiento.

- **Paquetes de R:**
  - `org.Hs.eg.db`
  - `DESeq2`
  - `ggplot2`
  - `clusterProfiler`
  - `dplyr`
  - `tidyr`
  - `readr`
  - Otros según necesidades del análisis

---


