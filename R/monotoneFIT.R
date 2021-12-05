monotoneFIT = function(mono.sim, grid.x = NULL, size.grid.x = 100){

  # If grid values of x are not provided, take a grid of size = size.grid.x with equispaced values from -1 and 1.
  if (is.null(grid.x)){
    grid.x = seq(-1, 1, length.out = size.grid.x)
  }
  else{
    # If provided, convert grid.x into a vector
    grid.x = as.vector(grid.x)
  }
  xi.mtx = as.matrix(mono.sim$xi)
  # return(xi.mtx)
  # Call c++ function g_mtx to calculate g(x) for the grid.x using xi's obtained from the MCMC.
  g.sample = g_mtx(xi_mtx = xi.mtx, grid_x = grid.x, u = mono.sim$knots)


}
