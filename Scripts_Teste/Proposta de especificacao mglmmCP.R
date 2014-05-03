## Testes com especificação do modelo -- Usuário especifica como um modelo 
## univariado e a função muda a fórmula para formato multivariado adequado.
# Edit 1: Talvez seja melhor pensar em uma forma de especificar em três partes:
#         1. Parâmetros que variam conforme a variável resposta;
#         2. Parâmetros com valor comum em todas as variáveis resposta;
#         3. Variável identificadora das variáveis resposta.


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