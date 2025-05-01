
#En este script (realizado en RStudio) se incluyen los comandos que han permitido el cálculo de los porcentajes de retención intrónica, la generación de boxplots y el análisis de rutas de enriquecimiento para
#los genes con mayor porcentaje de retención intrónica

#Paquetes requeridos para el análisis
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

#----------------------------------------------------------------------------------------------------------------
#CUANTIFICACIÓN: CÁLCULO DE PORCENTAJES DE RETENCIÓN INTRÓNICA
#----------------------------------------------------------------------------------------------------------------
#Introducir el input_dir donde están los archivos 
files <- list.files(input_dir, pattern = "_introns.txt", full.names = TRUE)

for (file in files) {
  sample_id <- gsub("_introns.txt$", "", basename(file))
  
  df <- read_delim(file, delim = " ", col_names = c("gene", "intron_count"), show_col_types = FALSE, trim_ws = TRUE)
  df <- mutate(df, ID = sample_id)
  
  dataframes[[sample_id]] <- df
  
}


# Combina todos los data.frames en uno solo
all_data_htseq <- bind_rows(dataframes)
View(all_data_htseq)

#Cálculo del porcentaje de retención intrónica: utiliza metadata (archivo SRARunTable)


all_data_htseqprueb <- left_join(all_data_htseq, metadata[, c("Run", "braak_lewy_body_stage")], by = c("ID" = "Run"))
View(all_data_htseqprueb)


summary_table <- all_data_htseqprueb %>%
  distinct(ID, braak_lewy_body_stage)

# Ver el resultado
View(summary_table)


all_data_htseqprueb2 <-all_data_htseqprueb %>%
  mutate(SampleID_Unique = paste0(braak_lewy_body_stage, "_", ID))
View(all_data_htseqprueb2)


yesno_htseq <- all_data_htseqprueb2 %>%
  group_by(SampleID_Unique) %>%
  summarise(no =sum(intron_count == 0),
            yes=sum(intron_count > 0),
            .groups = "drop") %>%
  mutate(percIR = yes / (no + yes) * 100)
View(yesno_htseq)

yesno_htseq_with_stage <- yesno_htseq %>%
  mutate(
    braak_lewy_body_stage = as.numeric(sub("^([0-9]+)_.*", "\\1", SampleID_Unique))  # Extrae el primer número
  )
View(yesno_htseq_with_stage)



percIR_by_stage <- yesno_htseq_with_stage %>%
  group_by(braak_lewy_body_stage) %>%
  summarise(
    percIR_avg = mean(percIR, na.rm = TRUE),  # Calcula el promedio de percIR por cada grupo
    .groups = "drop"  # Elimina el agrupamiento para un dataframe limpio
  )
View(percIR_by_stage)

#-------------------------------------------------------------------------------------------------------
#REPRESENTACIÓN DE BOXPLOTS PARA PORCENTAJE DE RETENCIÓN INTRÓNICA POR ETAPA DE BRAAK
#-------------------------------------------------------------------------------------------------------
#Si no asumimos normalidad:

#Kruskal-Wallis
kruskal_test <- kruskal.test(percIR ~ factor(braak_lewy_body_stage), data = yesno_htseq_with_stage)

# Ver resultados
print(kruskal_test)

#Como el p-value es menor que 0.05, hay diferencias significativas entre al menos un par de grupos
library(pgirmess)
kruskal_posthoc <- kruskalmc(yesno_htseq_with_stage$percIR, 
                             factor(yesno_htseq_with_stage$braak_lewy_body_stage))

# Ver los resultados
print(kruskal_posthoc)


# Cargar las librerías necesarias
library(ggplot2)
library(ggpubr)
library(pgirmess)

kruskal_posthoc <- kruskalmc(yesno_htseq_with_stage$percIR, factor(yesno_htseq_with_stage$braak_lewy_body_stage))

# Extraer las comparaciones significativas
signif_pairs <- kruskal_posthoc$dif.com[kruskal_posthoc$dif.com$stat.signif == TRUE, ]

# Crear una nueva columna para el número de asteriscos según la significancia
# Aquí, ya que stat.signif es TRUE, agregamos un asterisco por cada comparación significativa
signif_pairs$asterisks <- "*"

# Crear la lista de comparaciones para usar en ggplot
comparisons_list <- lapply(rownames(signif_pairs), function(x) strsplit(x, " - ")[[1]])

# Crear el gráfico con ggplot para los boxplots
ggplot(yesno_htseq_with_stage, aes(x = factor(braak_lewy_body_stage), y = percIR, fill = factor(braak_lewy_body_stage))) +
  geom_boxplot(outlier.shape = NA) +  # Crear los boxplots sin puntos atípicos
  geom_jitter(width = 0.2, alpha = 0.5) +  # Agregar los puntos para mejor visualización
  geom_signif(
    comparisons = list(c("0", "5")),  # Comparación entre los grupos 0 y 5
    annotations = "*",  # Coloca un asterisco para indicar significancia
    y_position = max(yesno_htseq_with_stage$percIR) + 5,  # Ajusta la posición vertical de la anotación
    map_signif_level = TRUE,  # Muestra el asterisco si es significativo
    textsize = 7,  # Tamaño del asterisco
    test = "wilcox.test"  # Método de prueba
  ) +
  labs(
       x = "Etapa de Braak ",
       y = "Porcentaje de Retención Intrónica (%)") +
  theme_minimal() +
  theme(legend.position = "none")  # Opcional: Eliminar la leyenda si no es necesaria


#---------------------------------------------------------------------------------------------------
#ANÁLISIS DE RUTAS DE ENRIQUECIMIENTO DE GENES
#--------------------------------------------------------------------------------------------------

# Obtener gene_symbol desde all_data_htseqprueb2
# Unimos los dataframes en la columna 'SampleID_Unique'

# Realizar la unión de los dataframes
yesno_htseq_with_stage_con_gen <- merge(yesno_htseq_with_stage, 
                                        all_data_htseqprueb2[, c("SampleID_Unique", "gene")], 
                                        by = "SampleID_Unique", 
                                        all.x = TRUE)
View(yesno_htseq_with_stage_con_gen)

# Aplicamos mapIds para obtener los símbolos de los genes
gene_symbol <- mapIds(
  org.Hs.eg.db,
  keys = yesno_htseq_with_stage_con_gen$gene,  # Utilizamos la columna 'gene' de tu dataframe
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"  # Tomamos el primer valor si hay múltiples símbolos asociados
)

# Añadimos la columna GeneSymbol a tu dataframe
yesno_htseq_with_stage_con_gen$GeneSymbol <- gene_symbol


# Obtener los primeros 150.000 genes ordenados por percIR
top150000_with_pathway <- yesno_htseq_with_stage_con_gen %>%
  filter(!is.na(GeneSymbol)) %>%
  arrange(desc(percIR)) %>%
  head(150000)


#Ordenamos por descendente
# Filtrar y ordenar usando dplyr
yesno_htseq_with_stage_con_gen <- yesno_htseq_with_stage_con_gen %>%
  filter(!is.na(GeneSymbol)) %>%  # Filtra filas donde GeneSymbol no es NA
  arrange(desc(percIR))  # Ordena por percIR en orden descendente


# Instalar clusterProfiler si no lo tienes y también org.Hs.eg.db
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")
BiocManager::install("clusterProfiler")

# Cargar el paquete
library(clusterProfiler)
# Cargar la librería para KEGG
library(org.Hs.eg.db)

# Convertir los nombres de los genes a identificadores ENSEMBL (o si tienes ENSEMBL, usa los nombres de genes directamente)
gene_list <- top150000_with_pathway$GeneSymbol

# Convertir los GeneSymbols a IDs de ENSEMBL usando org.Hs.eg.db
gene_entrez <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Filtrar las filas con un ENTREZID válido (número entero) de 4 dígitos
filtered_gene_entrez <- gene_entrez[nchar(as.character(gene_entrez$ENTREZID)) == 4, ]

# Verifica el resultado

head(filtered_gene_entrez)


# Realizar análisis de enriquecimiento para rutas de KEGG

kegg_enrich <- enrichKEGG(gene = filtered_gene_entrez$ENTREZID, organism = 'hsa')

# Ver los resultados del análisis
head(kegg_enrich)

