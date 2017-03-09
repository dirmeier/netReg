library(stringr)
library(ggplot2)
library(tidyr)
library(tibble)
library(microbenchmark)

setwd("~/PROJECTS/netreg_project/results/")

do.time <- function()
{
  bench.fls <- grep("time", list.files(".", full.names=T), value=T)
  bench.df <- c()
  ben.l <- list()
  for (fl in bench.fls)
  {
    m <- stringr::str_match(fl, "benchmark_time_n_([[:digit:]]+)_p_([[:digit:]]+)_.*")
    n <- as.numeric(m[2])
    p <- as.numeric(m[3])
    r <- readRDS(fl)
    bench.df <- rbind(bench.df, cbind(n,p,r))
    ben.l[[paste("n", n, "p", p, sep="_")]] <- r
  }
  df <- bench.df %>% dplyr::rename(N=n, P=p, Language=expr, Time=time) %>%
    dplyr::mutate(NP=paste(paste0("N=", N), paste0("P=", P), sep=", "),
                  N=as.factor(N),
                  P=as.factor(P),
                  Language=as.factor(Language))
  
  ben.l
  
}

do.rss <- function()
{
  bench.fls <- grep("rss", list.files(".", full.names=T), value=T)
  bench.df <- c()
  ben.l <- list()
  for (fl in bench.fls)
  {
    m <- stringr::str_match(fl, "benchmark_rss_n_([[:digit:]]+)_p_([[:digit:]]+)_q_([[:digit:]]+)_sig_([[:digit:]]+).*")
    n <- as.numeric(m[2])
    p <- as.numeric(m[3])
    q <- as.numeric(m[4])
    sig <- as.numeric(m[5])
    r <- readRDS(fl)
    bench.df <- rbind(bench.df, cbind(N=n, P=p ,Q=q, SD=sig, R=r))
    ben.l[[paste("n", n, "p", p, "q", q, "sig", sig, sep="_")]] <- r
  }
  df <- as_tibble(bench.df) %>% 
    tidyr::gather(Model, RSS, Lasso, Edgenet) %>%
    dplyr::mutate(Model=as.factor(Model),
                  MSE=log(RSS/N),
                  NPQ=as.factor(paste(paste0("N=", N), paste0("P=", P), paste0("Q=", Q), sep=", ")),
                  N=as.factor(paste0("N=", N)),
                  P=as.factor(paste0("P=", P)),
                  Q=as.factor(paste0("Q=", Q)),
                  sd=as.factor(SD))
  df$Noise <- "Low"
  df$Noise[df$SD == 2] <- "Medium"
  df$Noise[df$SD == 5] <- "High"
  df$Noise <- factor(df$Noise, levels=c("Low", "Medium", "High"))
  ggplot2::ggplot(df) +
    ggplot2::geom_boxplot(aes(x=Noise, y=MSE, fill=Model), width=0.5) +
    ggplot2::facet_grid(. ~ NPQ,  scales="free_y") +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.background=element_rect(fill="black")) +
    ggplot2::theme(strip.text=element_text(color="white", face="bold")) +
    ylab("log(MSE)")
}
