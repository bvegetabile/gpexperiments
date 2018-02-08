test_x <- matrix(seq(-3,3, length.out = 5000))
rbenchmark::benchmark(gpexperiments::sqexp_poly(test_x, c(1,1,1)),
                      par_sepkernel(test_x, c(1, 1), c(1, 1, 1)),
                      gpbalancer::par_sqexp(test_x, c(1,1,1)),
                      replications = 1)

all.equal(gpexperiments::sqexp_poly(test_x, c(1,1,4)),
          par_sepkernel(test_x, 
                        c(4, 1), 
                        c(1, 1, 1), sig_noise = 1e-4))

plot(gpexperiments::sqexp_poly(test_x, c(1,2,4)),
     par_sepkernel(test_x, c(4, 1), c(2, 1, 1)))

image(par_sepkernel(test_x, c(4, 1), c(2, 1, 1)))
