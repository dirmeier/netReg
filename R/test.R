library(tensorflow)

tfp <- reticulate::import("tensorflow_probability")


n <- 10000
x <- matrix(rnorm(n * 10), n)
b <- rnorm(10)
a <- 2
y <- a +  x %*% b + rnorm(n, 0, .1)

plot(x, y)

beta  <- tf$Variable(tf$zeros(shape(10L, 1L)))
alpha <- tf$Variable(tf$zeros(shape(1L)))

x <- tf$cast(x, tf$float32)
y <- tf$cast(y, tf$float32)

loss <- function(a, b) {
    mean <- tf$matmul(x, b) + a
    #di <- tfp$distributions$Normal(mean, st)
    #- di$log_prob(y) + + tf$abs(b)
    tf$reduce_mean(tf$pow(y - mean, 2)) + tf$abs(b) * 10
}

optimizer <- tf$train$GradientDescentOptimizer(learning_rate = 0.01)
train <- optimizer$minimize(loss(alpha, beta))

sess <- tf$Session()
sess$run(tf$global_variables_initializer())
for (step in seq(1000)) {
    sess$run(train)
}
cat(step, "->", sess$run(beta), sess$run(alpha), "\n")

