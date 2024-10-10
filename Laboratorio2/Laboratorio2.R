install.packages("tidyverse")
library(tidyverse)
data(starwars)



# Seleccionar todas las columnas menos el nombre
starwars %>% select(-name)
#Seleccionar sólo las columnas que tienen subraya (_)
starwars %>% select(contains("_"))
#Seleccionar sólo las columnas que empiezan con "s"
starwars %>% select(starts_with("s"))
#Crear un data frame con los nombres y planeta de origen (homeworld)
homeworld <- starwars %>% select(name, homeworld)
#Filtrar datos 
#Filtrar por especies: sólo humanos
human <- starwars %>% filter(species == "Human")
#Filtrar por especies: sólo humanos del planeta Tatooine
humantattoine <- starwars %>% filter(species == "Human", homeworld == "Tatooine")
#Crear un nuevo datframe con todas las especies menos los Droides
starwars_nodroids <- starwars %>% filter(species != "Droid")
#Usamos group_by y tally
n <- starwars %>% group_by(species)

#Añadiendo otra variable
starwars %>% group_by(species, gender) %>% tally()

#Si lo quieres guardar en el environment recuerda asignarle un nombre
table_gender <- starwars %>% group_by(species, gender) %>% tally()

species <- starwars %>% group_by(species) %>% summarise(mean_height = mean(height, na.rm = T),mean_mass = mean(mass,na.rm = T))





#EJERCICIO calcular la sd de esos parametros
speciessd <- starwars %>% group_by(species) %>% summarise(SD_height = sd(height, na.rm = T),sd_mass = sd(mass,na.rm = T)) 




#Hacer un gráfico de la altura vs. la masa de los personajes
ggplot(starwars, aes(height, mass)) + geom_point()

#Puedes modificar el color 
ggplot(starwars, aes(height, mass)) + geom_point(colour = "red")

#Modificando el color y el punto
ggplot(starwars, aes(height, mass)) + geom_point(colour = "purple", pch = 3)

#Modificando el color y el fondo 
ggplot(starwars, aes(height, mass)) + geom_point(colour = "red") + theme_light()






#EJERCICIO2
starwars %>% select("mass","name")
#jabba de hut pesa 1300kg
#quitamos a jaba y creamos el ggplot
starwars2 <- starwars %>% filter(name != "Jabba Desilijic Tiure")
ggplot(starwars2, aes(height, mass)) + geom_point(colour = "green") + theme_dark()




#EJERCICIO3

#CREA UN RESUMEN DE LA MEDIA
toy <- read_csv("Descargas/toy.csv")

Resumen <- toy %>% group_by(Sex) %>% summarise(mediaIMC = mean(IMC,na.rm = T),mediaIAS = mean(IAS,na.rm = T),mediapeso = mean(Weight_Kg,na.rm = T),mediaCcintura = mean(Ccintura,na.rm = T))


#SOLO PACIENTES FEMENINOS
femeninos <- toy %>% filter(Sex == "Women")
toy %>% filter(Sex == "Women") %>% tally
#58
femeninos %>% filter(IMC_clas == "Overweight") %>% tally()
#9

#Crea un gráfico IMC y Weight
ggplot(toy, aes(IMC,Weight_Kg)) + geom_point(colour = "green") + theme_dark()
obesidad <- toy %>% filter(IMC_clas != "Normal")
ggplot(obesidad, aes(IMC,Weight_Kg)) + geom_point(colour = "green") + theme_dark()





#Instalar los packages
install.packages("ape")
install.packages("phangorn")
install.packages("phytools")
