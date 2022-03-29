library(copula)
library(ggplot2)
Clayton <- claytonCopula(2)
n <- 2500
matUV.Clayton <- rCopula(n, Clayton)
dataUV.Clayton <- as.data.frame(matUV.Clayton)
colnames(dataUV.Clayton) <- c("u", "v")
plot(matUV.Clayton)
ggplot(data = dataUV.Clayton, aes(x = u, y = v)) + 
  geom_point(color = "slateblue", size = 2, alpha = 0.6) +
  xlab("u") + 
  ylab("v") +
  ggtitle("Clayton con teta = 2")+
  theme_minimal()
dev.new()
Gumbel <- gumbelCopula(2)
n <- 2500
matUV.Gumbel <- rCopula(n, Gumbel)
dataUV.Gumbel <- as.data.frame(matUV.Gumbel)
colnames(dataUV.Gumbel) <- c("u", "v")
plot(matUV.Gumbel)
ggplot(data = dataUV.Gumbel, aes(x = u, y = v)) + 
  geom_point(color = "slateblue", size = 2, alpha = 0.6) +
  xlab("u") + 
  ylab("v") +
  ggtitle("Gumbel con teta = 2")+
  theme_minimal()



Frank <- frankCopula(10)
n <- 2500
matUV.Frank <- rCopula(n, Frank)
dataUV.Frank <- as.data.frame(matUV.Frank)
colnames(dataUV.Frank) <- c("u", "v")
plot(matUV.Frank)
ggplot(data = dataUV.Frank, aes(x = u, y = v)) + 
  geom_point(color = "slateblue", size = 2, alpha = 0.6) +
  xlab("u") + 
  ylab("v") +
  ggtitle("Frank con teta = 10")+
  theme_minimal()



Frank2 <- frankCopula(-10)
n <- 2500
matUV.Frank2 <- rCopula(n, Frank2)
dataUV.Frank2 <- as.data.frame(matUV.Frank2)
colnames(dataUV.Frank2) <- c("u","v")
ggplot(data = dataUV.Frank2, aes(x = u, y = v)) + 
  geom_point(color = "slateblue", size = 2, alpha = 0.6) +
  xlab("u") + 
  ylab("v") +
  ggtitle("Frank con teta = -10")+
  theme_minimal()
plot(matUV.Frank2)