matrix <- read.delim("binary_formated.tsv")

matrix_new = t(matrix)

write.table(matrix_new, file = "matrix_transposed.tsv")
