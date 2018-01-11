n <- 2000
n_test <- 1000
x <- as.matrix(rnorm(n))
x_test <- as.matrix(seq(min(x), max(x), length.out = n_test))
y <- x^4 + 2*x + rnorm(n, sd = 5)

y_true <- x_test^4 + 2*x_test
out_ests <- gpr(x, y, sqexp, c(1), 1e-1, x_test)
plot(x, y, pch = 19, col=rgb(0,0,0,0.5))
points(x, out_ests$f_train, pch=4)
lines(x_test, out_ests$f_test)
lines(x_test, y_true, col='red')
rbenchmark::benchmark(gpr(x, y, sqexp, c(1), 1e-1, x_test), replications = 1)
plot(x_test, (y_true - out_ests$f_test))
abline(h=0)
