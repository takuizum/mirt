expect_class <- function(x, class) expect_true(inherits(x, class))

test_that('IRTpars', {

    dat <- expand.table(LSAT7)
    mod0 <- mirt(dat, 1, '3PL', verbose=FALSE)
    cfs <- coef(mod0)
    expect_equal(unname(cfs[[1L]][,'d']), 1.865620, tolerance = 1e-4)
    cfs <- coef(mod0, simplify = TRUE)
    expect_equal(unname(cfs$items[,'d']), c(1.86562047, 0.09404091, 1.78051762, 0.48189021, 1.86845282 ), tolerance = 1e-4)
    cfs <- coef(mod0, IRTpars = TRUE)
    expect_equal(unname(cfs[[1L]][,'b']), -1.852919, tolerance = 1e-4)
    cfs <- coef(mod0, IRTpars = TRUE, simplify = TRUE)
    expect_equal(unname(cfs$items[,'b']), c(-1.85291884, -0.04877029, -1.06840017, -0.65462572, -2.43551721), tolerance = 1e-4)

    mod1 <- mirt(dat, 1, '3PL', verbose=FALSE, SE=TRUE)
    cfs <- coef(mod1)
    expect_equal(unname(cfs[[1L]][,'d']), c(1.865620, 1.601845, 2.129396), tolerance = 1e-4)
    cfs <- coef(mod1, simplify = TRUE)
    expect_equal(unname(cfs$items[,'d']), c(1.86562047, 0.09404091, 1.78051762, 0.48189021, 1.86845282 ), tolerance = 1e-4)
    cfs <- coef(mod1, IRTpars = TRUE)
    expect_equal(unname(cfs[[1L]][,'b']), c(-1.852919, -2.362622, -1.343216), tolerance = 1e-4)
    cfs <- coef(mod1, IRTpars = TRUE, simplify = TRUE)
    expect_equal(unname(cfs$items[,'b']), c(-1.85291884, -0.04877029, -1.06840017, -0.65462572, -2.43551721), tolerance = 1e-4)
    cfs <- coef(mod1, printSE = TRUE)
    expect_equal(unname(cfs[[1L]][,'d']), c(1.8656205, 0.1345818), tolerance = 1e-4)
    expect_equal(unname(cfs[[2L]][,'logit(g)']), c(-0.8709602, 0.6723956), tolerance = 1e-4)
    cfs <- coef(mod1, IRTpars = TRUE, printSE = TRUE)
    expect_equal(unname(cfs[[1L]][,'b']), c(-1.8529188, 0.2600573), tolerance = 1e-4)
    expect_equal(unname(cfs[[2L]][,'g']), c(0.2950545, 0.1398565), tolerance = 1e-4)

    mod2 <- mirt(dat, 1, 'Rasch', verbose=FALSE)
    cfs <- coef(mod2, simplify=TRUE)
    expect_equal(unname(cfs$items[,'d']), c(1.8680155, 0.7908857, 1.4607798, 0.5214175, 1.9927099), tolerance = 1e-4)
    cfs <- coef(mod2, IRTpars=TRUE, simplify=TRUE)
    expect_equal(unname(cfs$items[,'b']), c(-1.8680155, -0.7908857, -1.4607798, -0.5214175, -1.9927099), tolerance = 1e-4)

    mod3 <- mirt(dat, 1, 'Rasch', verbose=FALSE, SE=TRUE)
    cfs <- coef(mod3)
    expect_equal(unname(cfs[[1L]][,'d']), c(1.868016, 1.671204, 2.064827), tolerance = 1e-4)
    cfs <- coef(mod3, printSE = TRUE)
    expect_equal(unname(cfs[[1L]][,'d']), c(1.868016, 0.100416), tolerance = 1e-4)
    cfs <- coef(mod3, IRTpars = TRUE)
    expect_equal(unname(cfs[[1L]][,'b']), c(-1.868072, -2.064887, -1.671256), tolerance = 1e-4)
    cfs <- coef(mod3, IRTpars = TRUE, printSE = TRUE)
    expect_equal(unname(cfs[[1L]][,'b']), c(-1.868016, 0.100416), tolerance = 1e-4)

})

test_that('IRTpars gpcm delta-method SEs', {

    # b_k = -(d_k - d_{k-1})/a with d_0 = 0 fixed; compare every b_k SE (b1 included)
    # against the delta method applied to vcov()
    mod <- mirt(Science, 1, 'gpcm', SE=TRUE, verbose=FALSE)
    cfs <- coef(mod, IRTpars = TRUE, printSE = TRUE)
    V <- vcov(mod)
    for(j in seq_len(ncol(Science))){
        it <- extract.item(mod, j)
        nms <- names(it@parnum)
        a <- it@par[nms == 'a1']
        ds <- it@par[grepl('^d[0-9]+$', nms)]
        lab <- paste0(nms, '.', it@parnum)
        pa <- lab[nms == 'a1']
        pd <- lab[grepl('^d[0-9]+$', nms)]
        for(k in 2L:length(ds)){
            if(k == 2L){
                grad <- c((ds[k] - ds[k-1L])/a^2, -1/a)
                ind <- c(pa, pd[k])
            } else {
                grad <- c((ds[k] - ds[k-1L])/a^2, 1/a, -1/a)
                ind <- c(pa, pd[k-1L], pd[k])
            }
            se <- sqrt(drop(t(grad) %*% V[ind, ind] %*% grad))
            expect_equal(unname(cfs[[j]]['SE', paste0('b', k-1L)]), unname(se), tolerance = 1e-6)
        }
    }

})
