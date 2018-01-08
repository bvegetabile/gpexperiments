
n <- 250
x <- as.matrix(rnorm(n))
x_test <- as.matrix(seq(min(x), max(x), length.out = 500))
y <- x^4 + 2*x + rnorm(n, sd = 5)

grid_res <- 25
sig_0 <- seq(1, 50, length.out = grid_res)
sig_1 <- seq(1, 50, length.out = grid_res)

ll_mat <- matrix(NA, grid_res, grid_res)
for(i in 1:grid_res){
    for(j in 1:grid_res){
        test_sig0 <- sig_0[i]
        test_sig1 <- sig_1[j]
        ll_mat[i, j] <- gpR(x, y, par_sqexp,
                            c(test_sig1, 5), test_sig0, x_test)[[2]]
    }
    message('.', appendLF = F)
}
image(sig_0, sig_1, ll_mat)
contour(sig_0, sig_1, ll_mat, add=T, nlevels = 5000)

max_ll <- which(ll_mat == max(ll_mat), arr.ind = T)
fin_sig0 <- sig_0[max_ll[1]]
fin_sig1 <- sig_1[max_ll[2]]

tester <- gpR(x, y, par_sqexp, c(fin_sig1, 2), fin_sig0, x_test)[[1]]

plot(x,y, pch=4)
points(x_test, tester, type='l')

