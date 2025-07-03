library(Matrix)
library(tidyr)
library(dplyr)
library(stringr)
library(purrr)
library(readr)
library(statmod)
library(tweedie)
library(corrplot)

load("RData/Rebuttal250121/CosMx1k6k_input_TweedieGLM.RData", verbose = T)

#### Identities
mData <- mDataIdent %>%
  mutate(Classif = ifelse(Classif == "OPC", "oligodendrocyte", Classif))
cellsByType <- cellsByTypeIdent

pData <- mData %>%
  right_join(cellsByType %>% as_tibble(rownames = "CellID"), by = "CellID")

pData <- pData %>%
  group_by(SampleID, FOV_ID, Platform, Classif) %>%
  mutate(nCells = n()) %>%
  relocate(nCells, .after = Type) %>%
  group_by(SampleID, FOV_ID, Platform, Classif, nCells) %>%
  summarise(across(GPMstructured:NEUdisp_core, ~ round(mean(.))), .groups = "drop")

pData <- pData %>%
  filter(nCells >= 5)

tmp <- mData %>%
  mutate(FOV_ID = factor(FOV_ID)) %>%
  filter(Type == "Malignant") %>%
  {table(.$FOV_ID, .$Classif)} %>%
  as.data.frame.matrix %>%
  .[, colnames(cellsByType)] %>%
  as_tibble(rownames = "FOV_ID")
colnames(tmp)[-1] <- str_c("n_", colnames(tmp)[-1])

pData <- pData %>%
  left_join(tmp, by = "FOV_ID")

pData <- pData %>%
  mutate(GPMstructured = ifelse(GPMstructured > n_GPMstructured, n_GPMstructured, GPMstructured),
         MTCstructured = ifelse(MTCstructured > n_MTCstructured, n_MTCstructured, MTCstructured),
         PPRstructured = ifelse(PPRstructured > n_PPRstructured, n_PPRstructured, PPRstructured),
         GPMdispersed = ifelse(GPMdispersed > n_GPMdispersed, n_GPMdispersed, GPMdispersed),
         MTCdispersed = ifelse(MTCdispersed > n_MTCdispersed, n_MTCdispersed, MTCdispersed),
         PPRdispersed = ifelse(PPRdispersed > n_PPRdispersed, n_PPRdispersed, PPRdispersed),
         NEUstructured = ifelse(NEUstructured > n_NEUstructured, n_NEUstructured, NEUstructured),
         NEUrim = ifelse(NEUrim > n_NEUrim, n_NEUrim, NEUrim),
         NEUdisp_core = ifelse(NEUdisp_core > n_NEUdisp_core, n_NEUdisp_core, NEUdisp_core)
  )

tmp <- pData[, 6:14]/pData[, 15:23]
colnames(tmp) <- str_c("r_", colnames(tmp))

pData <- pData %>%
  bind_cols(tmp)

pDataIdent <- pData
####

#### Activities
mData <- mDataAct %>%
  mutate(Classif = ifelse(Classif == "OPC", "oligodendrocyte", Classif))
cellsByType <- cellsByTypeAct

pData <- mData %>%
  right_join(cellsByType %>% as_tibble(rownames = "CellID"), by = "CellID")

pData <- pData %>%
  group_by(SampleID, FOV_ID, Platform, Classif) %>%
  mutate(nCells = n()) %>%
  relocate(nCells, .after = Type) %>%
  group_by(SampleID, FOV_ID, Platform, Classif, nCells) %>%
  summarise(across(GPMstructured:NEUdisp_core, ~ round(mean(.))), .groups = "drop")

pData <- pData %>%
  filter(nCells >= 5)

tmp <- mData %>%
  mutate(FOV_ID = factor(FOV_ID)) %>%
  filter(Type == "Malignant") %>%
  {table(.$FOV_ID, .$Classif)} %>%
  as.data.frame.matrix %>%
  .[, colnames(cellsByType)] %>%
  as_tibble(rownames = "FOV_ID")
colnames(tmp)[-1] <- str_c("n_", colnames(tmp)[-1])

pData <- pData %>%
  left_join(tmp, by = "FOV_ID")

pData <- pData %>%
  mutate(GPMstructured = ifelse(GPMstructured > n_GPMstructured, n_GPMstructured, GPMstructured),
         MTCstructured = ifelse(MTCstructured > n_MTCstructured, n_MTCstructured, MTCstructured),
         PPRstructured = ifelse(PPRstructured > n_PPRstructured, n_PPRstructured, PPRstructured),
         GPMdispersed = ifelse(GPMdispersed > n_GPMdispersed, n_GPMdispersed, GPMdispersed),
         MTCdispersed = ifelse(MTCdispersed > n_MTCdispersed, n_MTCdispersed, MTCdispersed),
         PPRdispersed = ifelse(PPRdispersed > n_PPRdispersed, n_PPRdispersed, PPRdispersed),
         NEUstructured = ifelse(NEUstructured > n_NEUstructured, n_NEUstructured, NEUstructured),
         NEUrim = ifelse(NEUrim > n_NEUrim, n_NEUrim, NEUrim),
         NEUdisp_core = ifelse(NEUdisp_core > n_NEUdisp_core, n_NEUdisp_core, NEUdisp_core)
  )

tmp <- pData[, 6:14]/pData[, 15:23]
colnames(tmp) <- str_c("r_", colnames(tmp))

pData <- pData %>%
  bind_cols(tmp)

pData <- pData %>%
  filter(Classif %in% c("Hypoxic_score", "Scavenger_Immunosuppressivescore", "Systemic_Inflammatoryscore", "Complement_Immunosuppressivescore", "Microglial_Inflammatoryscore"))

pDataAct <- pData
####

pData <- bind_rows(pDataIdent, pDataAct)

resList <- pData %>%
  split(.$Classif) %>%
  map(\(x){
    ddata <- x %>%
      dplyr::select(SampleID:NEUdisp_core) %>%
      pivot_longer(cols = 6:ncol(.), names_to = "Subtype", values_to = "n") %>%
      mutate(ID = str_c(FOV_ID, Subtype, sep = ";")) %>%
      relocate(ID, .before = SampleID)
    tmp <- x %>%
      dplyr::select(FOV_ID, contains("n_")) %>%
      pivot_longer(cols = 2:ncol(.), names_to = "Subtype", values_to = "N") %>%
      mutate(Subtype = gsub("n_", "", Subtype),
             ID = str_c(FOV_ID, Subtype, sep = ";")
      ) %>%
      dplyr::select(-FOV_ID, -Subtype)
    ddata <- ddata %>%
      left_join(tmp, by = "ID")
    tmp <- x %>%
      dplyr::select(FOV_ID, contains("r_")) %>%
      pivot_longer(cols = 2:ncol(.), names_to = "Subtype", values_to = "r") %>%
      mutate(Subtype = gsub("r_", "", Subtype),
             ID = str_c(FOV_ID, Subtype, sep = ";")
      ) %>%
      dplyr::select(-FOV_ID, -Subtype)
    ddata <- ddata %>%
      left_join(tmp, by = "ID")
    
    ddata <- ddata %>%
      filter(N > 0)
    
    
    (aGlm <- glm(n ~ Subtype - 1 + nCells, offset = log(N), family = tweedie(var.power = 1.9, link.power = 0), data = ddata))
    summary(aGlm)
    
    ccoef <- coef(aGlm)
    ccoef <- ccoef %>%
      {tibble(Subtype = gsub("Subtype", "", names(.)), betaValue = unname(.))} %>%
      mutate(Classif = ddata %>% slice(1) %>% pull(Classif)) %>%
      relocate(Classif, .before = Subtype)
    ccoef <- ccoef %>%
      slice(-grep("nCells", Subtype))
    
    tab <- tibble(Statistic = c("Deviance", "Pearson"),
                  GoF = c(aGlm$deviance, sum(aGlm$weights * aGlm$residuals^2)),
                  df = rep(aGlm$df.residual, 2),
                  pValue = NA
    )
    tab$pValue <- pchisq(tab$GoF, df = tab$df, lower.tail = FALSE)
    tab <- tab %>%
      mutate(Classif = ddata %>% slice(1) %>% pull(Classif)) %>%
      relocate(Classif, .before = Statistic)
    
    res <- list(ccoef = ccoef, tab = tab, ddata = ddata)
    return(res)
  })

ccoef <- resList %>%
  lapply(\(x) x$ccoef %>%
           mutate(betaValue = exp(betaValue)
           )
  ) %>%
  do.call(bind_rows, .)

ccoef <- ccoef %>%
  dplyr::select(Classif, Subtype, betaValue) %>%
  pivot_wider(names_from = Subtype, values_from = betaValue)
tmp <- ccoef %>%
  dplyr::select(-1) %>%
  as.matrix
rownames(tmp) <- ccoef$Classif
ccoef <- tmp

ccoef <- ccoef[c("Macrophage",  "Monocyte",  "Neutrophil", "cDC", "Microglia", "Tcells",
                 "Hypoxic_score",  "Scavenger_Immunosuppressivescore",  "Complement_Immunosuppressivescore", "Systemic_Inflammatoryscore", "Microglial_Inflammatoryscore",
                 "neurons", "astrocytes", "oligodendrocyte", "vasculature"),
               c("GPMstructured", "MTCstructured", "PPRstructured", "NEUstructured", "GPMdispersed", "MTCdispersed", "PPRdispersed", "NEUdisp_core", "NEUrim")]
# ccoef[is.na(ccoef)] <- 0

tab <- t(apply(ccoef, 1, \(x) (x - mean(x))/sd(x)))
tab[tab < -2] <- -2
tab[tab > 2] <- 2

rownames(tab) <- rownames(tab) %>%
  str_split("_") %>%
  map_chr(\(x) x[1])

corrplot(tab, method = "square", is.corr = F,
         insig = "blank",
         pch.cex = 1.8, pch.col = "white",
         col = rev(colorRampPalette(c("darkorange2", "tan2", "#FFD580", "white", "#4393C3", "#104A85", "#0D3D6E"))(51)),
         col.lim = c(-ceiling(max(tab)*10)/10, ceiling(max(tab)*10)/10),
         tl.col = "black", cl.pos = "b",  tl.srt = 35,  tl.cex = 1
)