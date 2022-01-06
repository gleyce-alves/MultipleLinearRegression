---
title: Regressão Linear Múltipla
author: "Gleyce Alves"
header-includes:
 - \usepackage{amsthm, amssymb}
 
output: html_document

fig_caption: true
---


```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```


<br/><br/>

<center> <h1> Regressão Linear Múltipla  no ![](http://developer.r-project.org/Logo/Rlogo-5.png) </h1> </center>

<br/><br/>

<div style='text-align: justify;'>

* O objetivo agora é mensurar o impacto sobre a média da variável resposta $y$ com mais de uma covariável.
O modelo linear de regressão múltipla
é dada por 
\begin{align*}
y_i &= 
\beta_1 + \beta_2x_{i2} + \beta_3x_{i3} + \cdots + \beta_kx_{iK} + \epsilon_i,
\quad i = 1, 2, \ldots, n.
\end{align*}


* Para exemplificar o modelo de $\textbf{Regressão Linear Múltipla}$
vamos utilizar o banco de dados $\textbf{prostate}$ disponível na biblioteca $\textit{faraway}$ do ambiente R. 
Os dados são de um estudo com 97 homens com câncer de próstata que
deveriam receber uma prostatectomia radical. Fornece uma comparação entre esses pacientes em termos do volume do câncer em $cm^3$ (lcavol),
logaritmo do peso da próstata em $g$ (lweight),
idade (age),
logaritmo da quantidade de hiperplasia prostática benigna (lbph),
invasão da vesícula seminal (1 = sim, 0 = não) (svi),
logaritmo da penetração capsular em $cm$ (lcp),
pontuação gleason (6, 7 ou > 7) (gleason), 
porcentagem da pontuação de Gleason 4 ou 5 (pgg45) e logaritmo do antígeno prostático específico em $ng/ml$ (lpsa).


+ O objetivo do modelo é estabelecer a relação entre a variável resposta $\textbf{lpsa}$ com as variáveis preditoras 
$\textbf{lcavol}$,
$\textbf{lweight}$, 
$\textbf{age}$,
$\textbf{lbph}$,
$\textbf{svi}$,
$\textbf{lcp}$,
$\textbf{gleason}$
e
$\textbf{pgg45}$.

$\textbf{Consideração:}$


Em geral, níveis mais altos de PSA são indicadores mais fortes de câncer de próstata. Um PSA de 4 corresponde a um lpsa de 1,39.

</div>

### -- Entrada dos dados -- Banco de dados
```{r}

library("faraway")
library("knitr")
dados <- prostate
kable(head(prostate), align = "c")
```

### -- Estatísticas descritivas
```{r}
summary(dados)
```

### -- Gráfico de dispersão
```{r}
plot(dados)
```

<div style='text-align: justify;'> 
$\textbf{Análise:}$
O gráfico de dispersão sugere que faz sentido modelar estas variáveis utilizando um modelo linear múltiplo. A regressão linear para este exemplo consiste basicamente de estimarmos os coeficientes $\beta_1, \beta_2,\dots, \beta_9$.
</div>

### -- Correlação entre as variáveis
```{r}
M <- cor(dados, method = "pearson")
library(corrplot)
library(RColorBrewer)
corrplot(M, type = "upper", method = "color", order = "original", outline = "white", 
        tl.col="black", tl.srt=45, addCoef.col = "black", tl.cex = 1.5, number.cex = 1.0,
        col = colorRampPalette(c("lightblue", "yellow", "brown"))(10))
```

<div style='text-align: justify;'> 
$\textbf{Análise:}$
Por meio do gráfico de correlação entre as variáveis explicativas é possível verificar que existem algumas variáveis relacionadas e que podem causar multicolineariedade. 
Veja por exemplo entre as variáveis $\textbf{lcavol}$ e $\textbf{svi}$.
Essa análise será feita mais adiante por meio do VIF (Fator de Inflação de Variância)
</div>


### -- Regressão Linear Múltipla
```{r}
# Auxiliares
nobs <- dim(dados)[1]
npar <- dim(dados)[2]

# matriz X
library("dplyr")
Covariaveis <- as.matrix(select(dados, -lpsa))
X <- cbind(1, Covariaveis)
Y <- dados[,9]

# vetor b
b <- solve(t(X)%*%X) %*% t(X)%*%Y

# Resíduos
echapeu <- Y - X%*%b
  
# Estimativa de sigma^2
sigma2chapeu <- crossprod(echapeu)/(nobs-npar)

# Erros-padrão
erroP <- sqrt(as.double(sigma2chapeu)*diag(solve(t(X)%*%X)))

# t-valor (Estatística de teste)
tvalor <- b/erroP

# p-valor
pvalor <- 2*(1-pt(abs(tvalor), df = nobs-npar))

# Medidas de Bondade de Ajuste
Yest  <- X%*%b
SQT   <- sum((Y-mean(Y))^2)
SQRes <- sum((Y - Yest)^2)
aux   <- (nobs - 1)/(nobs - npar)

# Coeficiente de determinação
R2   <- 1 - SQRes/SQT
R2a  <- 1 - aux*(1-R2) ## R2 ajustado
R2nc <- 1 - crossprod(erroP)/crossprod(Y) ## R2 não ajustado

# Critérios de informações
AIC <- nobs*log(SQRes/nobs) + 2*npar
BIC <- nobs*log(SQRes/nobs) + npar*log(nobs)/nobs

# Impressão dos dados
imp <- cbind(b, erroP, tvalor, pvalor)
colnames(imp) <- c("Estimativas", "DesvioP", "t-valor", "p-valor")
rownames(imp) <- c("Intercepto", "lcavol", "lweight", "age", "lbph", "svi", "lcp", "gleason", "pgg45")

Cba <- rbind(R2, R2a, R2nc, AIC, BIC)
colnames(Cba) <- c("Critérios de seleção")
rownames(Cba) <- c("R2", "R2a", "R2nc", "AIC", "BIC")

kable(round(imp, 4), align = "c")
kable(round(Cba, 4), align = "c")
```

<div style='text-align: justify;'> 
$\textbf{Análise:}$
O resultado do teste de hipóteses acima mostra a estimativa de cada coeficiente associada ao respectivo teste de hipótese. Como em outros testes de hipóteses, rejeitamos $\mathcal{H}_0$ se o valor de $p$ associado ao coeficiente for igual o menor que o nível de significância adotado. 
Neste exemplo, o valor de $p$ associado aos $\beta$'s foi maior que 0,05 (5\%)
com exceção das varáveis volume do câncer, logaritmo do peso da próstata e invasão da vesícula seminal.
Diante disso, não rejeitamos a hipótese nula,
ou seja, os $\beta$'s são iguais a zero com exceção dos $\beta$'s associados as variveis $\textbf{lcavol}$, $\textbf{lweight}$ e $\textbf{svi}$. Note que os coeficientes de determinação estão dando valores razoavelmente bons. Vamos verificar se há multicolinearidade entre as variáveis.
</div>


<left> <h3> 
No $\texttt{R}$ os modelos lineares são resolvidos com a função $\texttt{lm}$.
</h3> </left>


### -- Estimar os parâmetros do modelo -- M1
```{r}
M1 <- lm(lpsa ~ lcavol + lweight + age + lbph + svi + lcp + gleason + pgg45, data = dados)
M1a <- summary(M1)
summary(M1)
```


### -- Multicolinearidade
```{r}
library("regclass", quietly = TRUE)
kable(VIF(M1), align = 'c', col.names = "VIF")
```


<div style='text-align: justify;'> 
$\textbf{Análise:}$
Observamos que o VIF de todas as variáveis é menor que 5, ou seja, não há grandes evidências de ter multicolinearidade. 
</div>

### -- Estimar os parâmetros do modelo -- M2
```{r}
M2 <- lm(lpsa ~ lcavol + lweight + age + lbph + svi + lcp, data = dados)
M2a <- summary(M2)
summary(M2)
```

<div style='text-align: justify;'> 
$\textbf{Análise:}$
O resultado do teste de hipóteses acima mostra as estimativas de cada coeficiente associada ao respectivo teste de hipótese. Como em outros testes de hipóteses, rejeitamos $\mathcal{H}_0$ se o valor de $p$ associado ao coeficiente for igual o menor que o nível de significância adotado. 
Neste modelo, sem as variáveis $\textbf{gleason}$ e $\textbf{pgg45}$, os valores de $p$ associados a $\beta_2$, $\beta_3$ e $\beta_6$ foram menores que 0,05 (5\%) o que é uma forte evidência de que o volume do câncer, o logarítmo do peso da próstata e a invasão da vesícula seminal, de fato, são significantes. 
</div>

### -- Multicolinearidade
```{r}
kable(VIF(M2), align = 'c', col.names = "VIF")
```

<div style='text-align: justify;'> 
$\textbf{Análise:}$
Observamos que o VIF de todas as variáveis é menor que 5, ou seja, não há grandes evidências de ter multicolinearidade. 
</div>

### -- Estimar os parâmetros do modelo -- M3
```{r}
M3 <- lm(lpsa ~ lcavol + lweight + lbph + svi, data = dados)
M3a <- summary(M3)
summary(M3)
```

<div style='text-align: justify;'> 
$\textbf{Análise:}$
Nesse terceiro modelo retiramos as variáveis $\textbf{age}$, $\textbf{lcp}$, $\textbf{gleason}$ e $\textbf{pgg45}$. Assim, observamos que às variáveis 
$\textbf{lcavol}$, $\textbf{lweight}$ e $\textbf{svi}$ são significantes.
</div>

### -- Multicolinearidade
```{r}
kable(VIF(M3), align = 'c', col.names = "VIF")
```

<div style='text-align: justify;'> 
$\textbf{Análise:}$
Observamos que o VIF de todas as variáveis é menor que 5, ou seja, não há grandes evidências de ter multicolinearidade. 
</div>

### -- Estimar os parâmetros do modelo -- M4
```{r}
M4 <- lm(lpsa ~ lcavol + lweight + svi, data = dados)
M4a <- summary(M4)
summary(M4)
```

<div style='text-align: justify;'> 
$\textbf{Análise:}$
Nessa quarta modelagem deixamos apenas as variáveis 
$\textbf{lcavol}$, $\textbf{lweight}$ e $\textbf{svi}$
pelo fato de serem significantes, a 95% de confiança.
A seguir ilustramos o VIF dessas variáveis,
na qual foram abaixo de 5,
ou seja,
há grandes evidências de não ter multicolinearidade.
</div>

### -- Multicolinearidade
```{r}
kable(VIF(M4), align = 'c', col.names = "VIF")
```


### -- Comparação dos modelos
```{r}
## R2 não centrado
R2n1 <- 1 - crossprod(M1$resid)/crossprod(Y)
R2n2 <- 1 - crossprod(M2$resid)/crossprod(Y)
R2n3 <- 1 - crossprod(M3$resid)/crossprod(Y)
R2n4 <- 1 - crossprod(M4$resid)/crossprod(Y)
  
Mod1 <- c(AIC(M1), BIC(M1), M1a$r.squared, M1a$adj.r.squared, R2n1)
Mod2 <- c(AIC(M2), BIC(M2), M2a$r.squared, M2a$adj.r.squared, R2n2)
Mod3 <- c(AIC(M3), BIC(M3), M3a$r.squared, M3a$adj.r.squared, R2n3)
Mod4 <- c(AIC(M4), BIC(M4), M4a$r.squared, M4a$adj.r.squared, R2n4)

aux <- rbind(Mod1, Mod2, Mod3, Mod4)
colnames(aux) <- c("AIC", "BIC", "R2", "R2 Ajustado", "R2 Não centrado")
kable(aux, align = "c", digits = 4)
```

### -- Teste Reset
#### -- Para saber se precisa fazer transformação nos dados, ou seja, para verificar se o modelo está correto  ou não na especificação usamos o teste $\texttt{Reset}$. As hipoteses são:

\begin{align*}
\mathcal{H}_0 &: \text{o modelo corretamente especificado}
\quad \text{versus} \quad
\mathcal{H}_1 : \text{o modelo incorretamente especificado}
\end{align*}

```{r}
try(require(lmtest, quietly = TRUE), silent=TRUE)
pvalor <- cbind(resettest(M1, type = "regressor")[4], resettest(M2, type = "regressor")[4],
resettest(M3, type = "regressor")[4], resettest(M4, type = "regressor")[4])

colnames(pvalor) <- c("Modelo 1", "Modelo 2", "Modelo 3", "Modelo 4")
rownames(pvalor) <- "p-valor"

kable(pvalor, align = "c")
```


### -- Teste da razão de verossimilhança
```{r}
Ma1 <- lm(lpsa ~ svi, data = dados)
Ma2 <- lm(lpsa ~ lcavol + lweight + svi, data = dados)
Ma3 <- lm(lpsa ~ lcavol + lweight + lbph + svi, data = dados)
Ma4 <- lm(lpsa ~ lcavol + lweight + age + lbph + svi + lcp, data = dados)
Ma5 <- lm(lpsa ~ lcavol + lweight + age + lbph + svi + lcp + gleason + pgg45, data = dados)
anova(Ma1, Ma2, Ma3, Ma4, Ma5)
```


<div style='text-align: justify;'> 
$\textbf{Análise:}$
Após a estimação dos modelos e verificado se as covariáveis existem multicolinearidades,
fizemos a comparação dos 4 modelos através dos coeficientes de determinações
e critérios de informações alternativos. Além disso,
aplicamos o teste Reset de Ramsey e o teste da razão de verossimilhança. Pelo teste da razão de verossimilhança não foi possível concluir qual o modelo melhor especificado. Por outro lado, pelo teste Reset e de acordo com os valores de AIC, BIC, R2, R2 Ajustado e R2 Não centrado, obtidos na comparação dos modelos, é possível concluir que o modelo 4 se mostrou melhor, mais eficiente.
Agora vamos fazer o estudo de diagnóstico do modelo final (modelo 4).
Mas, antes imprimiremos as estimativas e algumas ferramentas básicas do Modelo Final, o quarto modelo.
</div>

<center> <h3> Modelo Final </h3> </center>

### -- Estimativas
```{r}
MF <- lm(lpsa ~ lcavol + lweight + svi, data = dados)
summary(MF)
```
```{r}
kable(anova(MF), align = "c", digits = 4)
```

### -- 

<div style='text-align: justify;'> 
$\textbf{Análise:}$
Como o p-valor é menor que $\alpha$ (considerando $\alpha = 0,05$), rejeitamos a hipótese nula e com isso temos que o modelo é bem explicado, ou seja,
as variáveis preditoras está bem correlacionada com a variável resposta. 

O coeficiente de determinação ajustado foi de 0,6144, indicando que 61,44% da variabilidade dos dados é explicada pelo modelo. Ressalta-se que para modelos de regressão que contêm diferentes números de preditores, utilizamos o coeficiente de determinação ajustado. 


Sendo assim temos o seguinte modelo:
\begin{align*}
\text{lpsa} &= -0,26809 + (0,55164)\text{lcavol} + (0,50854)\text{lweight} + (0,66616)\text{svi}.
\end{align*}

Quando o peso da próstata e a invasão da vesícula seminal são mantidas constantes, 0,55164 é uma estimativa do aumento esperado do antígeno prostático específico (PSA), quando há um aumento no volume do câncer.

Agora quando temos o volume do câncer e o peso da próstata constantes, 0,66616 uma estimativa do aumento esperado do antígeno prostático específico (PSA), quando há um aumento na invasão da vesícula seminal.



Por último, consideramos que não há invasão de vesícula seminal (svi=0), e consideramos os valores médios de lweight e lcavol, obtemos lpsa sendo 2,3725.
</div>



<br/><br/>

<center> <h3> Diagnóstico </h3> </center>


-- Quando terminamos de ajustar o modelo de regressão devemos avaliar se estes pressupostos são válidos. No modelo de regressão normal linear, assumimos independência dos erros, homocedasticidade (variância constante dos erros) e que os erros têm distribuição normal.

-- Para verificar essas suposições existem diversos métodos gráficos e testes estatísticos. Aqui vamos apresentar algumas estatísticas e gráficos que ajudam na identificação de possíveis problemas na regressão.


### -- Teste Reset
```{r}
resettest(MF, type = "regressor")
```

### -- Multicolinearidade
```{r}
kable(VIF(MF), align = 'c', col.names = "VIF")
```
### -- Normalidade do resíduos

Primeiro, faremos o teste de $\textbf{Shapiro Wilk}$ para analisar a normalidade dos resíduos:

$H_0$: distribuição dos dados $=$ normal $\rightarrow p>0,05$

$H_1$: distribuição dos dados $\neq$ normal $\rightarrow p \leq 0,05$



```{r}
shapiro.test(MF$residuals)
```

Nesse caso, como o $p$ foi maior que 0,05 podemos consideram a normalidade dos dados.

### -- Teste de Bera-Jarque para Normalidade

Temos que:

$H_0$: distribuição dos dados $=$ normal $\rightarrow p>0,05$

$H_1$: distribuição dos dados $\neq$ normal $\rightarrow p \leq 0,05$

```{r}
library(normtest, quietly = TRUE)
library(MASS, quietly = TRUE)
sresid <- studres(MF) 
ajb.norm.test(sresid)
```
Pelo Teste de Bera-Jarque também não rejeitamos a hipótese nula, logo, consideramos a normalidade dos dados.

### -- Teste da razão de verossimilhança
```{r}
MFa  <- lm(lpsa ~ lcavol + lweight + svi, data = dados)
anova(MFa)
```


<div style='text-align: justify;'> 
</div>




<br/><br/>



<justify> <h4> 
-- Com o foco no diagnóstico, será verificado a normalidade dos resíduos studentizados por meio do qqplot e histograma.
</h3> </justify>




### -- qqPlot 
```{r}
qqnormb <- function(X, Y, Z){
li <- qqnorm(X, plot.it = F)
x  <- qqnorm(Y, plot.it = F)
ls <- qqnorm(Z, plot.it = F)
 
data <- data.frame(lix = li$x, liy = li$y, xx= x$x, xy = x$y ,lsx = ls$x, lsy = ls$y)
 
require(ggplot2)
 p <- ggplot(data, aes(lix,liy, xx,xy, lsx,lsy)) + geom_ribbon(aes(ymin = liy, ymax = lsy), alpha = 0.3, col = "grey") + geom_point(aes(xx, xy))
 return(p)
}
#---------------------------------------------------------------#
par(mfrow=c(1,1))
X <- model.matrix(MF)
n <- nrow(X)
p <- ncol(X)
H <- X%*%solve(t(X)%*%X)%*%t(X)
h <- diag(H)
si <- lm.influence(MF)$sigma
r <- resid(MF)
tsi <- r/(si*sqrt(1-h))
#
ident <- diag(n)
epsilon <- matrix(0,n,100)
e <- matrix(0,n,100)
e1 <- numeric(n)
e2 <- numeric(n)
#
for(i in 1:100){
 epsilon[,i] <- rnorm(n,0,1)
 e[,i] <- (ident - H)%*%epsilon[,i]
 u <- diag(ident - H)
 e[,i] <- e[,i]/sqrt(u)
 e[,i] <- sort(e[,i]) }
#
for(i in 1:n){
 eo <- sort(e[i,])
 e1[i] <- (eo[2]+eo[3])/2
 e2[i] <- (eo[97]+eo[98])/2 }
#
med <- apply(e,1,mean)
faixa <- range(tsi,e1,e2)
#
#--------------------------------------------------------------#
 
qqnormb(e1, tsi, e2) + theme_bw() + xlab("Percentil da Normal(0,1)") + ylab("Resíduo Studentizado") + geom_abline(slope = 1, intercept = 0, col = "blue")
```

<div style='text-align: justify;'> 
$\textbf{Análise:}$
O gráfico qqplot dos resíduos com banda de confiança serve para checar a normalidade dos erros. Se todos os pontos estiverem dentro da banda, podemos assumir que os erros têm distribuição normal.  Sendo assim verificamos que todos os pontos estão sobre a reta e dentro das bandas com 95% de confiança. A seguir temos o histograma da distribuição dos resíduos Studentizados para confirmar essa hipótese.
</div>


### - Histograma dos Resíduos Studentizados
```{r}
library(MASS)
sresid <- studres(MF) 
hist(sresid, freq=FALSE, main="Distribuição dos Resíduos Studentizados", xlab = "Resíduos Studentizados", ylab = "Densidade", xlim = c(-2.8, 2.8))
xm <- seq(-1+min(sresid), 1+max(sresid), length=40) 
ym <- dnorm(xm) 
lines(xm, ym)
```




<br/><br/>




### -- Medidas de Alavanca, Pontos Influentes, Discrepantes e Outliers
#### -- Comandos auxiliares para os demais gráficos
```{r}
X <- model.matrix(MF)  # Matriz do modelo ajustado
lms <- summary(MF)
n <- nrow(X)
p <- ncol(X)
H <- X%*%solve(t(X)%*%X)%*%t(X) # Matriz Hat
h <- diag(H)  # Diagonal da matriz Hat (hii, medida de alavanca)
s <- lms$sigma
r <- resid(lms) # Residuos
ts <- r/(s*sqrt(1-h))
di <- (1/p)*(h/(1-h))*(ts^2) # Distância de cook
si <- lm.influence(MF)$sigma
tsi <- r/(si*sqrt(1-h)) # Resíduo Studentizado
```


#### -- Agora vamos ao primeiro gráfico, o da medida de alavanca
```{r}
h <- data.frame(h = h, i = c(1:(length(h))))
ggplot(h, aes(i, h)) + geom_point() + ylim(0, 0.4) + xlab("Índice") + ylab("Medida h") + theme_bw() + geom_abline(slope = 0, intercept = 2*p/n, col = "blue", type = "dashed") + geom_abline(slope = 0, intercept = 3*p/n, col = "red", type = "dashed") + annotate("text", x=62, y=0.26, label= "2k/n", col = "blue") + annotate("text", x=62, y=0.39, label= "3k/n", col = "red")
```

<div style='text-align: justify;'> 
$\textbf{Análise:}$ Esse gráfico serve para identificar pontos de alavanca, que podem influenciar desproporcionalmente na inclinação da reta ajustada.
Neste caso percebemos que os pontos, exceto um, estão abaixo da reta mostrando está no padrão.
</div>


<br/><br/>

#### -- Agora vamos fazer o gráfico da Distância de Cook,
```{r}
di <- data.frame(di = di, i = c(1:97))
plot <- ggplot(di, aes(i, di)) + geom_point() + xlab("Índice") + ylab("Distância de Cook") + theme_bw() + geom_abline(slope = 0, intercept = 8/(nobs - 2*4), col = "blue")
plot + annotate("text", x=22, y=0.32, label= "8/(n-2k)", col = "blue")
```

<div style='text-align: justify;'> 
$\textbf{Análise:}$ Esse gráfico ajuda a identificar pontos influentes, que causam alterações substanciais nos parâmetros quando são retirados do ajuste.
Neste caso percebemos que os pontos estão abaixo da reta mostrando está no padrão.
</div>




<br/><br/>


#### -- Agora vamos fazer o gráfico de resíduos studentizados pelo índice
```{r}
tsi <- data.frame(tsi = tsi, i = 1:nobs)
ggplot(tsi, aes(i, tsi)) + geom_point() + theme_bw() + xlab("Índice") + ylab("Resíduo Studentizados") + geom_hline(yintercept = -2, col = "blue") + geom_hline(yintercept = 2, col = "blue")
```


<div style='text-align: justify;'> 
$\textbf{Análise:}$ O gráfico de resíduos studentizados pelo índice serve para verificar se a suposição de independência dos erros é plausível. Além disso serve para identificar pontos discrepantes, com resíduo alto.
Neste caso percebemos que os pontos estão entre $-2$ e 2 com exceção de cinco pontos fora do intervalo $[-2,2]$. Pórem, mesmo assim, se mostra no padrão.
</div>




<br/><br/>

#### -- O próximo gráfico é de resíduo studentizado versus valor ajustado

```{r}
df <- data.frame(tsi = tsi$tsi, fit = MF$fitted)
ggplot(df, aes(fit, tsi)) + geom_point() + geom_hline(yintercept = -2, col = "blue") + geom_hline(yintercept = 2, col = "blue") + xlab("Valores Ajustados") + ylab("Resíduo Studentizado") + theme_bw()
```


<div style='text-align: justify;'> 
$\textbf{Análise:}$ O gráfico de resíduos studentizados pelos valores ajustados serve para verificar se a suposição de independência dos erros é plausível. Além disso serve para identificar pontos discrepantes, com resíduo alto.
Neste caso percebemos que os pontos estão entre $-2$ e 2 com exceção de 5 pontos fora do intervalo $[-2,2]$. Mas, mesmo assim, se mostra no padrão.
</div>




<br/><br/>

### -- DFFITS
```{r}
library(stats)
dff <- data.frame(dff = dffits(MF), i = 1:nobs)
plot <- ggplot(dff, aes(i, dff)) + geom_point() + theme_bw() + xlab("Índice") + ylab("DFFITS") + geom_hline(yintercept = 2*sqrt(4/nobs), col = "blue")
plot + annotate("text", x=25, y=0.80, label= expression(2*sqrt(4/n)), col = "blue")
```

<div style='text-align: justify;'> 
$\textbf{Análise:}$ O gráfico DFFITS serve para mostrar como é influente o ponto em uma regressão linear múltipla.
Neste caso percebemos que os pontos estão abaixo da reta com exceção de 5 pontos que estão acima da reta. Mas, mesmo assim, consideramos no padrão.
</div>



<br/><br/>


