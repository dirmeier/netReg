library(stringr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)
library(xtable)
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

  id <- group_indices(df, N, P, Language)
  df$grp <- id
  tab <- do.call("rbind", lapply(unique(id), function(e) {
        g <- dplyr::filter(df, grp==e)
        times <- quantile(g$Time, probs=c(0.25, 0.5, 0.75))
        cbind(N=as.character(g$N[1]), P=as.character(g$P[1]),
              Language=as.character(g$Language)[1],
              Min=min(g$Time),Mean=mean(g$Time), Median=median(g$Time), Max=max(g$Time))
    })) %>% as.data.table %>%
      dplyr::transmute(n=as.integer(N),
                    p=as.integer(P),
                    Model=as.character(Language),
                    Min=as.double(Min) / 1000 / 1000,
                    Max=as.double(Max) / 1000 / 1000,
                    Median=as.double(Median)   / 1000 / 1000,
                    Mean=as.double(Mean)   / 1000 / 1000)
  align <- rep("r", ncol(tab) + 1)
  bold <- function(x) {paste('{\\textbf{',x,'}}', sep ='')}
  caption <- "Timings of pure R vs R/C++ implementation in seconds."
  label <- "tab:timings"
  print(xtable(tab, align=align, digits=2,
               label=label, caption=caption),
        comment=F, booktabs=T,
        sanitize.colnames.function=bold,
        sanitize.subheadings.function = bold,
        include.rownames = FALSE)

df}

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
                  MS=(RSS/N),
                  NPQ=as.factor(paste(paste0("n=", N), paste0("p=", P), paste0("q=", Q), sep=", ")),
                  n=as.factor(paste0("N=", N)),
                  p=as.factor(paste0("P=", P)),
                  q=as.factor(paste0("Q=", Q)),
                  sd=as.factor(SD))
  df$Noise <- "Low"
  df$Noise[df$SD == 2] <- "Medium"
  df$Noise[df$SD == 5] <- "High"
  df$Model <- as.character(df$Model)
  df$Model[df$Model == "Edgenet"] <- "NetReg"
  df$Noise <- factor(df$Noise, levels=c("Low", "Medium", "High"))
  g <- ggplot2::ggplot(df) +
    ggplot2::geom_boxplot(aes(x=Noise, y=MSE, fill=Model), width=0.5) +
    ggplot2::facet_grid(. ~ NPQ,  scales="free_y") +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.background=element_rect(fill="black")) +
    ggplot2::theme(strip.text=element_text(color="white", face="bold", size=15)) +
    ggplot2::theme(text=element_text(size=16),
                   axis.text.x = element_text(size=16, angle=25),
                   legend.text=element_text(size=15)) +
    ylab("Mean Residual Sum of Squares")
  ggsave(g, filename="~/Desktop/plot.pdf", width=8, height=5)
}
