# lab weight
E7 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7_112018_Eim_COMPLETE.csv")
E10 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E10_112020_Eim_COMPLETE.csv")
CLS <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/CLS_complete.csv")
E11 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E11_012021_Eim_record.csv")

# select weight stuff

E7w <- select(E7, EH_ID, labels, experiment, dpi, relative_weight, mouse_strain, primary_infection, challenge_infection, 
              infection_history, OPG, batch)
E10w <- select(E10, EH_ID, labels, experiment, dpi, relative_weight, mouse_strain, primary_infection, challenge_infection, 
               infection_history, OPG, batch)
CLSw <- select(CLS, EH_ID, labels, experiment, dpi, relative_weight, mouse_strain, primary_infection, challenge_infection, 
               infection_history, OPG, batch)
E11w <- select(E11, EH_ID, labels, experiment, dpi, relative_weight, mouse_strain, primary_infection, challenge_infection, 
               infection_history, OPG, batch)
# add empty OPG for E11 until wehave it
E11w$OPG <- NA

W <- rbind(E7w, E10w)
W <- rbind(W, CLSw)
W <- rbind(W, E11w)
# fix mouse strains function


W$mouse_strain[W$mouse_strain == "SCHUNT"] <- "SCHUNT_SCHUNT"
W$mouse_strain[W$mouse_strain == "PWD"] <- "PWD_PWD"
# W$mouse_strain[W$mouse_strain == "BUSNA_PWD"] <- "PWD_BUSNA"
# W$mouse_strain[W$mouse_strain == "STRA_BUSNA"] <- "BUSNA_STRA"
# W$mouse_strain[W$mouse_strain == "PWD_SCHUNT"] <- "SCHUNT_PWD"
# W$mouse_strain[W$mouse_strain == "STRA_SCHUNT"] <- "SCHUNT_STRA"
W$simple_strain[W$mouse_strain == "SWISS"] <- "SWISS"
W$simple_strain[W$mouse_strain == "SCHUNT_SCHUNT"] <- "Mmd"
W$simple_strain[W$mouse_strain == "PWD_PWD"] <- "Mmm"
W$simple_strain[W$mouse_strain == "BUSNA_STRA"] <- "Hybrid"
W$simple_strain[W$mouse_strain == "STRA_BUSNA"] <- "Hybrid"
W$simple_strain[W$mouse_strain == "PWD_SCHUNT"] <- "Hybrid"
W$simple_strain[W$mouse_strain == "STRA_STRA"] <- "Mmd"
W$simple_strain[W$mouse_strain == "STRA_SCHUNT"] <- "Mmd"
W$simple_strain[W$mouse_strain == "PWD_BUSNA"] <- "Mmm"
W$simple_strain[W$mouse_strain == "SCHUNT_PWD"] <- "Hybrid"
W$simple_strain[W$mouse_strain == "SCHUNT_STRA"] <- "Mmd"
W$simple_strain[W$mouse_strain == "BUSNA_BUSNA"] <- "Mmm"
W$simple_strain[W$mouse_strain == "BUSNA_PWD"] <- "Mmm"

W$eimeria_species_primary[W$primary_infection == "E64"] <- "FER"
W$eimeria_species_challenge[W$challenge_infection == "E64"] <- "FER"
W$eimeria_species_primary[W$primary_infection == "E88"] <- "FAL"
W$eimeria_species_challenge[W$challenge_infection == "E88"] <- "FAL"
W$eimeria_species_primary[W$primary_infection == "Eflab"] <- "FAL"
W$eimeria_species_challenge[W$challenge_infection == "Eflab"] <- "FAL"
W$eimeria_species_primary[W$primary_infection == "E139"] <- "FER"
W$eimeria_species_challenge[W$challenge_infection == "E139"] <- "FER"
W$eimeria_species_primary[W$primary_infection == "UNI"] <- "UNI"
W$eimeria_species_challenge[W$challenge_infection == "UNI"] <- "UNI"

W$eimeria_species_history[W$infection_history == "UNI:UNI"] <- "UNI:UNI"
W$eimeria_species_history[W$infection_history == "E88:E88"] <- "FAL:FAL"
W$eimeria_species_history[W$infection_history == "E88:E64"] <- "FAL:FER"
W$eimeria_species_history[W$infection_history == "E88:UNI"] <- "FAL:UNI"
W$eimeria_species_history[W$infection_history == "Eflab:UNI"] <- "FAL:UNI"
W$eimeria_species_history[W$infection_history == "Eflab:E88"] <- "FAL:FAL"
W$eimeria_species_history[W$infection_history == "E64:E64"] <- "FER:FER"
W$eimeria_species_history[W$infection_history == "E64:E88"] <- "FER:FAL"
W$eimeria_species_history[W$infection_history == "E64:UNI"] <- "FER:UNI"
W$eimeria_species_history[W$infection_history == "UNI:E88"] <- "UNI:FAL"
W$eimeria_species_history[W$infection_history == "UNI:E64"] <- "UNI:FER"
W$eimeria_species_history[W$infection_history == "E139:E64"] <- "FER:FER"
W$eimeria_species_history[W$infection_history == "E139:E88"] <- "FER:FAL"
W$eimeria_species_history[W$infection_history == "E139:UNI"] <- "FER:UNI"
W$eimeria_species_history[W$infection_history == "Eflab:E64"] <- "FAL:FER"


W$relevant_history[W$eimeria_species_history == "UNI:UNI"] <- "UNI"
W$relevant_history[W$eimeria_species_history == "UNI:FAL"] <- "FAL"
W$relevant_history[W$eimeria_species_history == "FAL:UNI"] <- "FAL"
W$relevant_history[W$eimeria_species_history == "FER:UNI"] <- "FER"
W$relevant_history[W$eimeria_species_history == "UNI:FER"] <- "FER"
W$relevant_history[W$eimeria_species_history == "FAL:FAL"] <- "FAL:FAL"
W$relevant_history[W$eimeria_species_history == "FER:FER"] <- "FER:FER"
W$relevant_history[W$eimeria_species_history == "FER:FAL"] <- "FER:FAL"
W$relevant_history[W$eimeria_species_history == "FAL:FER"] <- "FAL:FER"

write.csv(W, "C:/Users/exemp/OneDrive/Documents/W.csv")
# primary weightloss
ggplot(subset(W, !is.na(W$primary_infection)), aes(x = dpi, relative_weight, color = eimeria_species_primary)) +
  geom_jitter(width = 0.2) +
  geom_smooth() +
  facet_wrap(~simple_strain) +
  ggtitle(label = "Weight loss in primary infection")
# primaryoocyst shedding
ggplot(subset(W, !is.na(W$primary_infection)), aes(x = dpi, OPG, color = eimeria_species_primary)) +
  geom_jitter(width = 0.2) +
  geom_smooth() +
  facet_wrap(~simple_strain) +
  ggtitle(label = "oocyst shedding in primary infection")
# challenge weight loss
ggplot(subset(W, !is.na(W$challenge_infection)), aes(x = dpi, relative_weight, color = eimeria_species_challenge)) +
  geom_jitter(width = 0.2) +
  geom_smooth() +
  facet_wrap(~simple_strain) +
  ggtitle(label = "Weight loss in challenge infection")


ggplot(subset(W, !is.na(W$challenge_infection)), aes(x = dpi, relative_weight, color = eimeria_species_history)) +
  geom_jitter(width = 0.2) +
  geom_smooth() +
  facet_grid(eimeria_species_history~simple_strain) +
  ggtitle(label = "Weight loss in challenge infection with history")

# challenge oocysts
ggplot(subset(W, !is.na(W$challenge_infection)), aes(x = dpi, OPG, color = eimeria_species_challenge)) +
  geom_jitter(width = 0.2) +
  geom_smooth() +
  facet_wrap(~simple_strain) +
  ggtitle(label = "oocyst shedding in challenge infection")

ggplot(subset(W, !is.na(W$challenge_infection)), aes(x = dpi, OPG, color = eimeria_species_history)) +
  geom_jitter(width = 0.2) +
  geom_smooth() +
  facet_grid(eimeria_species_history~simple_strain) +
  ggtitle(label = "oocyst shedding in challenge infection with history")

# correlate weight at dissection vs maximum weight loss
swap <- read.csv("C:/Users/exemp/OneDrive/Documents/pca.data.swap.csv")

xx <- swap$relative_weight
yy <- swap$maximum_weight_loss_challenge

cor(xx, yy, method = "pearson", use = "complete.obs")

ggplot(swap, aes(x = relative_weight, y = maximum_weight_loss_challenge)) +
  geom_point()
