library(dplyr)
library(ggplot2)
library(gridExtra)

dados <- read.csv("dados.csv", sep=";")

chikv <- subset(dados, chikv_igm %in% c("Positivo", "Negativo"))
table(chikv$chikv_igm, chikv$sexo)

chikv_data <- data.frame(
  chikv_igm = c("Negativo", "Negativo", "Positivo", "Positivo"),
  sexo = c("Feminino", "Masculino", "Feminino", "Masculino"),
  count = c(45, 94, 13, 21),
  percent = c(0.33 , 0.67, 0.38, 0.62)
)

grafico <- ggplot(chikv_data, aes(x = chikv_igm, y = sexo, fill = percent)) +
  geom_tile(color = "black", size = 1) +
  geom_text(aes(label = percent), color = "black", size = 5, fontface = "bold") +  
  scale_fill_gradient(low = "#f8f8fa", high = "#7f0000") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(face = "bold", size = 8),
        axis.title = element_text(face = "bold", size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("CHIKV IgM") +  
  ylab("Sexo") +
  labs(fill = "Percentual")
#-------------------------------------------------------------------------------

denv <- subset(dados, den_igm %in% c("Positivo", "Negativo"))
table(denv$den_igm, denv$sexo)

denv_data <- data.frame(
  den_igm = c("Negativo", "Negativo", "Positivo", "Positivo"),
  sexo = c("Feminino", "Masculino", "Feminino", "Masculino"),
  count = c(38, 104, 8, 10),
  percent = c(0.27 , 0.73, 0.45, 0.55)
)

grafico2 <- ggplot(denv_data, aes(x = den_igm, y = sexo, fill = percent)) +
  geom_tile(color = "black", size = 1) +
  geom_text(aes(label = percent), color = "black", size = 5, fontface = "bold") +  
  scale_fill_gradient(low = "#f8f8fa", high = "#7f0000") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(face = "bold", size = 8),
        axis.title = element_text(face = "bold", size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("DENV IgM") +  
  ylab("Sexo") +
  labs(fill = "Percentual")

grafico3 <- grid.arrange(grafico, grafico2, ncol = 2)
ggsave("heatmap_genero.jpg", plot = grafico3, width = 12, height = 4, dpi = 600)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
table(chikv$chikv_igm, chikv$faixa_etaria)

chikv_idade <- data.frame(
  chikv_igm = c("Negativo", "Negativo", "Positivo", "Positivo"),
  idade = c("Jovem adulto (18-35 anos)", "Adulto (36-60 anos)", 
           "Jovem adulto (18-35 anos)", "Adulto (36-60 anos)"),
  count = c(78, 61, 25, 8),
  percent = c(0.56 , 0.44, 0.76, 0.24)
)

grafico4 <- ggplot(chikv_idade, aes(x = chikv_igm, y = idade, fill = percent)) +
  geom_tile(color = "black", size = 1) +
  geom_text(aes(label = percent), color = "black", size = 5, fontface = "bold") +  
  scale_fill_gradient(low = "#f8f8fa", high = "#7f0000") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(face = "bold", size = 8),
        axis.title = element_text(face = "bold", size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("CHIKV IgM") +  
  ylab("Faixa etária") +
  labs(fill = "Percentual")

#-------------------------------------------------------------------------------
table(denv$den_igm, denv$faixa_etaria)

denv_idade <- data.frame(
  denv_igm = c("Negativo", "Negativo", "Positivo", "Positivo"),
  idade = c("Jovem adulto (18-35 anos)", "Adulto (36-60 anos)", 
            "Jovem adulto (18-35 anos)", "Adulto (36-60 anos)"),
  count = c(84, 57, 8, 10),
  percent = c(0.60 , 0.40, 0.45, 0.55)
)

grafico5 <- ggplot(denv_idade, aes(x = denv_igm, y = idade, fill = percent)) +
  geom_tile(color = "black", size = 1) +
  geom_text(aes(label = percent), color = "black", size = 5, fontface = "bold") +  
  scale_fill_gradient(low = "#f8f8fa", high = "#7f0000") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(face = "bold", size = 8),
        axis.title = element_text(face = "bold", size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("DENV IgM") +  
  ylab("Faixa etária") +
  labs(fill = "Percentual")

grafico6 <- grid.arrange(grafico4, grafico5, ncol = 2)
ggsave("heatmap_idade.jpg", plot = grafico6, width = 14, height = 4, dpi = 600)