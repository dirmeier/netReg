library(stringr)
library(ggplot2)
library(microbenchmark)

setwd("~/PROJECTS/netreg_project/results/")
bench.fls <- list.files(".", full.names=T)
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
bench.df

ben.l

bench.df$expr <- factor(bench.df$expr)

ggplot2::ggplot(bench.df) +
  ggplot2::geom_boxplot(aes(as.factor(p), time, colour=as.factor(expr))) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)))  +
  ggplot2::theme_bw() 