library(readxl)

xl <- read_excel('Results/Heiser_results_5.0.xlsx', col_names=F, skip=3, sheet=1)
xl <- xl[,1:2]
xl <- xl[apply(is.na(xl), 1, sum)<2,]
xl$X0 <- sub('H2O ', '', xl$X0)
xl$X0 <- sub('Neural net', 'N.N.', xl$X0)
xl$X0 <- sub('RF', 'R.F.', xl$X0)
xl$X0 <- sub('BaTFLED3D ', '', xl$X0)
xl$X0 <- sub('Cl mean', 'Cl. mean', xl$X0)
xl$X0 <- sub('Dr mean', 'Dr. mean', xl$X0)
xl$X0 <- sub('Cl dr mean', 'Cl. Dr. mean', xl$X0)
write.table(xl, file='Results/CV_design_matrix.txt', 
            quote=F, row.names=F, sep='\t')

xl <- read_excel('Results/Heiser_results_5.0.xlsx', col_names=F, skip=3, sheet=2)
xl <- xl[,1:2]
xl <- xl[apply(is.na(xl), 1, sum)<2,]
xl$X0 <- sub('H2O ', '', xl$X0)
xl$X0 <- sub('Neural net', 'N.N.', xl$X0)
xl$X0 <- sub('RF', 'R.F.', xl$X0)
xl$X0 <- sub('BaTFLED3D ', '', xl$X0)
xl$X0 <- sub('Cl mean', 'Cl. mean', xl$X0)
xl$X0 <- sub('Dr mean', 'Dr. mean', xl$X0)
xl$X0 <- sub('Cl dr mean', 'Cl. Dr. mean', xl$X0)
write.table(xl, file='Results/CV_kern_design_matrix.txt', 
            quote=F, row.names=F, sep='\t')

xl <- read_excel('Results/Heiser_results_5.0.xlsx', col_names=F, skip=3, sheet=3)
xl <- xl[,1:2]
xl <- xl[apply(is.na(xl), 1, sum) < 2,]
xl$X0 <- sub('H2O ', '', xl$X0)
xl$X0 <- sub('Neural net', 'N.N.', xl$X0)
xl$X0 <- sub('RF', 'R.F.', xl$X0)
xl$X0 <- sub('BaTFLED3D ', '', xl$X0)
xl$X0 <- sub('Cl mean', 'Cl. mean', xl$X0)
xl$X0 <- sub('Dr mean', 'Dr. mean', xl$X0)
xl$X0 <- sub('Cl dr mean', 'Cl. Dr. mean', xl$X0)
write.table(xl, file='Results/Final_design_matrix.txt', 
            quote=F, row.names=F, sep='\t')

xl <- read_excel('Results/Heiser_results_5.0.xlsx', col_names=F, skip=3, sheet=4)
xl <- xl[,1:2]
xl <- xl[apply(is.na(xl), 1, sum) < 2,]
xl$X0 <- sub('H2O ', '', xl$X0)
xl$X0 <- sub('Neural net', 'N.N.', xl$X0)
xl$X0 <- sub('RF', 'R.F.', xl$X0)
xl$X0 <- sub('BaTFLED3D ', '', xl$X0)
xl$X0 <- sub('Cl mean', 'Cl. mean', xl$X0)
xl$X0 <- sub('Dr mean', 'Dr. mean', xl$X0)
xl$X0 <- sub('Cl dr mean', 'Cl. Dr. mean', xl$X0)
write.table(xl, file='Results/Final_kern_design_matrix.txt', 
            quote=F, row.names=F, sep='\t')
