yeojohnson <- function(mat,
                       prop = 0.7,
                       ...){
    requireNamespace("recipes")
    requireNamespace("rsample") 
    split <- rsample::initial_split(data.frame(as.matrix(mat)),
                                    prop = prop)
    train  <- rsample::training(split)
    test   <- rsample::testing(split)
    rec <- recipes::recipe(train)
    yj_transform <- recipes::step_YeoJohnson(recipe = rec,
                                             recipes::all_numeric(),
                                             ...)
    yj_estimates <- recipes::prep(yj_transform, training = train)
    norm_dat <- function(dat,
                         yj_estimates){
        yj_te <- recipes::bake(yj_estimates, dat)
        mat_yj <- Matrix::Matrix(as.matrix(yj_te), sparse = TRUE)
        rownames(mat_yj) <- rownames(dat)
        return(mat_yj)
    }
    mat_yj_tr <- norm_dat(dat = train, yj_estimates = yj_estimates)
    mat_yj_te <- norm_dat(dat = test, yj_estimates = yj_estimates)
    mat_yj <- rbind(mat_yj_tr, mat_yj_te)
    return(mat_yj)
}
