#### Variables

# NK var
df <- read.csv("data/NKgenes_list.csv", sep=',')

# Create new columns for synonyms gene names
df_expanded <- df %>%
  separate(synonyms, into = paste0("synonyms_", 1:13), sep = ",", 
           fill = "right", extra = "drop")

# Remove white spaces at the beginning and at the end of each column value
df_expanded <- as.data.frame(lapply(df_expanded, str_trim), stringsAsFactors = FALSE)
df_expanded <- df_expanded %>% relocate(X, .before = 3)

df_expanded <- df_expanded %>%
  mutate(geneset = case_when(
    grepl("Inh", geneset) ~ gsub("Inh", "Inhibitory", geneset),
    grepl("Act", geneset) ~ gsub("Act", "Activating", geneset),
    TRUE ~ geneset
  ))

res_df <- data.frame(geneset = df_expanded$geneset, 
                     Hugo_Symb = df_expanded$HUGO_SYMBOL, 
                     final_symb = rep(0, nrow(df_expanded)),
                     is_synonym_used = rep(NA, nrow(df_expanded)))





