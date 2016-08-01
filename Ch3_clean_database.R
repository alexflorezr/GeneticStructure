str(DATABASE)
DATABASE_bckup <- DATABASE
ch3_database <- DATABASE
ch3_database <- ch3_database[-which(ch3_database$Species == "Alces_alces"),]
ch3_database <- ch3_database[-which(ch3_database$Species == "Balaena_mysticetus"),]
ch3_database <- ch3_database[-which(ch3_database$Species == "Canis_lupus"),]
ch3_database <- ch3_database[-which(ch3_database$Species == "Capra_pyrenaica"),]
ch3_database <- ch3_database[-which(ch3_database$Species == "Castor_fiber"),]
ch3_database <- ch3_database[-which(ch3_database$Species == "Crocuta_crocuta"),]
ch3_database <- ch3_database[-which(is.na(ch3_database$Latitude)),]
ch3_database <- ch3_database[-which(ch3_database$Latitude == ""),]
ch3_database <- ch3_database[-which(ch3_database$Mean_Age > 50000),]


ch3_database[which(ch3_database$Mean_Age > 50000)]
