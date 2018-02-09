test_x <- matrix(seq(-3,3, length.out = 500))
rbenchmark::benchmark(gpexperiments::sqexp_poly(test_x, c(1,1,1)),
                      par_sepkernel(test_x, c(1, 1, 1, 1, 1)),
                      gpbalancer::par_sqexp(test_x, c(1,1,1)),
                      replications = 10)

all.equal(gpexperiments::sqexp_poly(test_x, c(1,1,1)),
          par_sepkernel(test_x,
                        c(1, 1, 1, 1, 1), sig_noise = 1e-4))

plot(gpexperiments::sqexp_poly(test_x, c(1,2,4)),
     par_sepkernel(test_x, c(4, 1, 2, 1, 1)))

image(test_x, test_x, par_sepkernel(test_x, c(.25, 1, 2, 1, 1)), col = heat.colors(1000))
contour(test_x, test_x, par_sepkernel(test_x, c(.25, 1, 2, 1, 1)), add=T)
