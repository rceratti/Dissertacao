## Testes com especifica��o do modelo -- Usu�rio especifica como um modelo 
## univariado e a fun��o muda a f�rmula para formato multivariado adequado.
# Edit 1: Talvez seja melhor pensar em uma forma de especificar em tr�s partes:
#         1. Par�metros que variam conforme a vari�vel resposta;
#         2. Par�metros com valor comum em todas as vari�veis resposta;
#         3. Vari�vel identificadora das vari�veis resposta.


value~-1+variable+period:variable+(-1+variable|ID)

    main= value~period+(1|id)
multiVar=      ~variable
  common=       NULL


value~-1+variable+period+(-1+variable|ID)

    main= value~(1|id)
multiVar=      ~variable
  common=      ~period


til.1<-strsplit(paste(main),'\\~')




form<-resp~trt+tempo*trt+(1|id1)+(tes|id2)
id<-'composto'

form.ch<-paste(form)

pred.ch<-form.ch[grep('\\+',form.ch)]
pred.ch.1<-strsplit(pred.ch,"\\+")[[1]]

re.ind<-grep("\\|",pred.ch.1)
fix.ch<-pred.ch.1[-re.ind]
re.ch<-pred.ch.1[re.ind]

paste(pred.ch.1,":",id)

til<-form.ch[grep('\\~',form.ch)]

resp~-1+trt:composto+tempo:composto+(-1+composto|id)