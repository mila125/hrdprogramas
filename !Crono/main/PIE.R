CEP_csv  <- read.csv2("C:\\Users\\milax\\Documents\\!Crono_November_2020\\Data base\\ActualDatabase.csv")

for (j in 1:i) 
{
  aux<-scan("", what = character(),1);   #leo los  elementos que deben ser incluidos en el compuesto
  if(!(CEP_csv %>% dplyr::filter(Name == aux)))  #Averiguo sieste elemento existe
  {
    printf("try again!");
    j<-j-1; #si el elemento no existe en la db , el Crono lo pide ingressar nuevamente
  }
}

