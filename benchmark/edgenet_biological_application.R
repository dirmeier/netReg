library(data.table)
library(igraph)
library(magrittr)
library(netReg)
library(ggplot2)
library(stringr)
library(dplyr)
library(datastructures)

options(stringsAsFactors=FALSE)

.mse <- function(X.text, Y.test, coefs, intercepts)
{
    pr <- X.test %*% coefs + rep(1, nrow(X.test)) %*% t(intercepts)
    mse <- mean((Y.test - pr)**2)
    mse
}

fls        <- list.files("/Users/simondi/PROJECTS/netreg_project/results/", pattern="out_", full.names=T)
data.files <- "/Users/simondi/PROJECTS/netreg_project/code/netReg/data/"

fls <- c(list.files("/Users/simondi/PROJECTS/netreg_project/results/", full.names=T,
                    pattern="out_lasso_lambda_10_fold_"),
         list.files("/Users/simondi/PROJECTS/netreg_project/results/", full.names=T,
                    pattern="out_netReg_lambda_10_xi_4_")
)

tab <- NULL
for (fl in fls)
{

    mat <- stringr::str_match(fl, "out_(\\w+)_(lambda)_(\\d+.?\\d?)_(?:(xi)_(\\d+.?\\d*)_)?fold_(.+)_(\\w+).tsv")
    model <- stringr::str_to_title(mat[2])
    lambda <- as.numeric(mat[4])
    xi <- as.numeric(mat[6])
    tab <- rbind(tab,
                 data.frame(Model=model,
                            Lambda=lambda,
                            Xi=xi,
                            Fold=as.integer(mat[7]),
                            Parameter=mat[8],
                            Filename=fl))
}

tab

tab <- as.data.table(tab) %>%
    distinct(Model, Lambda, Xi, Fold, Parameter, .keep_all=TRUE)
grps <- dplyr::group_indices(tab, Model, Lambda, Xi)
tab$GRP <- grps


MSE <- NULL
for (g in grps)
{
    for (i in seq(10))
    {
        grp.dat <- dplyr::filter(tab, GRP == g, Fold==i)
        coefs   <- as.matrix(read.csv(as.character(dplyr::filter(grp.dat, Parameter=="coefficients")$Filename), sep=",", header=F))
        intercepts   <- as.matrix(read.csv(as.character(dplyr::filter(grp.dat, Parameter=="intercepts")$Filename), sep=",", header=F))
        X.test  <- as.matrix(read.csv(paste0(data.files, "X_test_", i, ".tsv"), sep="\t", header=F))
        Y.test  <- as.matrix(read.csv(paste0(data.files, "Y_test_", i, ".tsv"), sep="\t", header=F))
        mse <- .mse(X.test, Y.test, coefs, intercepts)
        MSE <- rbindlist(list(MSE,
                     data.table(Model=grp.dat$Model[1],
                                Fold=i,
                                Lambda=grp.dat$Lambda[1],
                                Xi=grp.dat$Xi[1],
                                MSE=mse)))
    }
    print(g)
}

MSE <- saveRDS(MSE, "./benchmark/mse_best_two.rds")
MSE <- readRDS( "./benchmark/mse_best_two.rds")

summ.mse <-
    group_by(MSE, Model, Lambda, Xi) %>%
    dplyr::summarize(mse=mean(MSE)) %>%
    dplyr::group_by(Model) %>%
    dplyr::slice(which.min(mse))

summ.mse

lam  <- filter(MSE, (Lambda==10 & Model=="Lasso")) %>% distinct(Model, Lambda, Fold, .keep_all=T)
grap <- filter(MSE, Lambda==10 , Xi==4) %>% distinct(Model, Lambda, Fold, .keep_all=T)
dat <- rbind(lam, grap)

dat$Model[dat$Model == "Netreg"] <- "NetReg"

g <-
    ggplot2::ggplot(dat) +
    ggplot2::geom_boxplot(aes(x=Model, y=MSE, fill=Model), width=0.5) +
    ggplot2::theme_bw() +
    coord_flip() +
    ggplot2::theme(strip.background=element_rect(fill="black")) +
    ggplot2::theme(strip.text=element_text(color="white", face="bold", size=15)) +
    ggplot2::scale_x_discrete(limits=c("NetReg", "Lasso")) +
    ggplot2::theme(text=element_text(size=16),
                   axis.text.x = element_text(size=16),
                   axis.text.y = element_text(size=16),
                   legend.text=element_text(size=17)) +
    xlab("") +
    ylab("Mean Residual Sum of Squares")

for (d in c("~/PROJECTS/netreg_project/plots/biological_application.eps",
            "~/PROJECTS/netreg_project/docs/paper/resubmission/biological_application.eps"))
{
    ggsave(g, filename=d,
           width=8, height=4, device="eps")
}
