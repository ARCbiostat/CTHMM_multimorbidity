
library(MMLCA)
library(ggplot2)

snack_nhm$Age_disc<- snack_nhm$Age/3
snack_nhm$MP_label <- factor(as.factor(snack_nhm$MP), levels = as.character(c(1,2,3)), labels = c("Mild MM", "Complex MM","Death"))

all_snack <- ggalluvial(
  snack_nhm,
  time_var = "Age_disc",
  "lopnr",
  "MP_label",
  colors= c("lightblue", "dodgerblue", "grey"),
  space = 20,
  "Death"
)

all_snack$plot+scale_x_continuous("Age",breaks = seq(20,35,by=2),labels = c("60-63","66-69", "72-75", "78-81", "84-87", "90-93", "96-99", "100+"))
